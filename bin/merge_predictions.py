#!/usr/bin/env python3
"""Merge TD2, MetaEuk, and GeneMarkS-T ORF predictions into a single best-per-gene set.

Usage:
    merge_predictions.py --td2-faa X --td2-map X \
        --metaeuk-faa X --metaeuk-map X --metaeuk-psauron X \
        --gmst-faa X --gmst-map X --gmst-psauron X \
        --salmon-quants dir1,dir2 \
        --min-psauron 0.3 --species SPECIES

Merge strategy (per gene):
    1. Drop genes with zero expression across all samples.
    2. Collect predictions from all sources that pass PSAURON >= min_psauron.
       Rescue: if no prediction passes but the gene has mean TPM >= 1,
       accept the best prediction even below the threshold (expression rescue).
    3. Rank by: completeness → PSAURON → expression (mean TPM) → protein length.
    4. Pick the top-ranked prediction.

Outputs:
    - {species}.faa       : merged protein FASTA (one per gene)
    - merge_map.tsv       : per-gene merge details
    - merge_stats.txt     : summary statistics for thinning report
"""

import argparse
import csv
import os
import sys
from Bio import SeqIO


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def parse_psauron_csv(path):
    """Parse a PSAURON output CSV (handles preamble lines).

    PSAURON uses different column names depending on input mode:
      - Protein mode (-p):    ``in-frame_score``  (hyphen)
      - Nucleotide CDS mode:  ``in_frame_score``  (underscore)
    We accept both.
    """
    scores = {}
    with open(path) as f:
        for line in f:
            if line.startswith("description,"):
                break
        else:
            print(f"WARNING: no 'description,' header in {path}", file=sys.stderr)
            return scores
        columns = line.strip().split(",")
        reader = csv.DictReader(f, fieldnames=columns)
        # Detect score column (hyphen vs underscore)
        if "in-frame_score" in columns:
            score_col = "in-frame_score"
        elif "in_frame_score" in columns:
            score_col = "in_frame_score"
        else:
            print(f"WARNING: no score column found in {path} (columns: {columns})",
                  file=sys.stderr)
            return scores
        for row in reader:
            try:
                scores[row["description"]] = float(row[score_col])
            except (ValueError, KeyError):
                continue
    return scores


def classify_completeness(seq):
    """Heuristic completeness from protein sequence (MetaEuk / fallback)."""
    if seq.startswith("M"):
        return "complete"
    return "partial"


COMPL_RANK = {"complete": 2, "3prime_partial": 1, "5prime_partial": 1,
              "partial": 0, "internal": 0, "unknown": 1}


def median(vals):
    if not vals:
        return 0
    s = sorted(vals)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def mean(vals):
    return sum(vals) / len(vals) if vals else 0


RESCUE_TPM_THRESHOLD = 1.0   # mean TPM to rescue sub-threshold genes


def parse_quant_sf(path):
    """Parse a Salmon quant.sf file → dict {Name: TPM}."""
    tpms = {}
    with open(path) as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.split("\t")
            tpms[parts[0]] = float(parts[3])   # Name, Length, EffLen, TPM, NumReads
    return tpms


def load_expression(quant_dirs):
    """Load Salmon quant.sf from each directory → per-gene mean TPM.

    Parameters
    ----------
    quant_dirs : str
        Comma-separated list of directories, each containing ``quant.sf``.

    Returns
    -------
    dict  {gene_id: mean_tpm}
    """
    all_tpms = {}          # gene_id → [tpm_sample1, tpm_sample2, ...]
    dirs = [d.strip() for d in quant_dirs.split(",") if d.strip()]
    n_samples = 0
    for d in dirs:
        sf = os.path.join(d, "quant.sf")
        if not os.path.isfile(sf):
            print(f"WARNING: {sf} not found, skipping", file=sys.stderr)
            continue
        n_samples += 1
        for gene_id, tpm in parse_quant_sf(sf).items():
            all_tpms.setdefault(gene_id, []).append(tpm)
    # Genes missing from a sample get TPM=0 for that sample
    gene_means = {}
    for gene_id, tpm_list in all_tpms.items():
        # pad with zeros for samples where the gene was absent
        while len(tpm_list) < n_samples:
            tpm_list.append(0.0)
        gene_means[gene_id] = sum(tpm_list) / n_samples
    print(f"Expression: loaded {n_samples} samples, {len(gene_means):,} genes with TPM data")
    return gene_means


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="3-way ORF prediction merge")
    ap.add_argument("--td2-faa", required=True)
    ap.add_argument("--td2-map", required=True)
    ap.add_argument("--metaeuk-faa", required=True)
    ap.add_argument("--metaeuk-map", required=True)
    ap.add_argument("--metaeuk-psauron", required=True)
    ap.add_argument("--gmst-faa", required=True)
    ap.add_argument("--gmst-map", required=True)
    ap.add_argument("--gmst-psauron", required=True)
    ap.add_argument("--salmon-quants", default=None,
                    help="Comma-separated Salmon quant directories")
    ap.add_argument("--min-psauron", type=float, default=0.3)
    ap.add_argument("--species", required=True)
    args = ap.parse_args()

    min_ps = args.min_psauron
    species = args.species

    # --- Load expression data (if provided) ---
    gene_tpms = {}
    has_expression = False
    if args.salmon_quants:
        gene_tpms = load_expression(args.salmon_quants)
        has_expression = bool(gene_tpms)

    # --- Parse TD2 proteins and scores ---
    td2_seqs = {}
    for rec in SeqIO.parse(args.td2_faa, "fasta"):
        td2_seqs[rec.id] = str(rec.seq)

    td2_scores = {}
    td2_compl = {}
    with open(args.td2_map) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            gene_id = parts[0]
            td2_scores[gene_id] = float(parts[2])
            # completeness column added by updated select_best_orf.py
            td2_compl[gene_id] = parts[4] if len(parts) > 4 else "unknown"

    print(f"TD2: {len(td2_seqs):,} proteins loaded")

    # --- Parse MetaEuk proteins and PSAURON scores ---
    metaeuk_seqs = {}
    for rec in SeqIO.parse(args.metaeuk_faa, "fasta"):
        metaeuk_seqs[rec.id] = str(rec.seq)

    metaeuk_scores = parse_psauron_csv(args.metaeuk_psauron)
    print(f"MetaEuk: {len(metaeuk_seqs):,} proteins, {len(metaeuk_scores):,} PSAURON scores")

    # --- Parse GeneMarkS-T proteins and scores ---
    gmst_seqs = {}
    for rec in SeqIO.parse(args.gmst_faa, "fasta"):
        gmst_seqs[rec.id] = str(rec.seq)

    gmst_scores = parse_psauron_csv(args.gmst_psauron)

    gmst_compl = {}
    with open(args.gmst_map) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                gmst_compl[parts[0]] = parts[2]

    print(f"GeneMarkS-T: {len(gmst_seqs):,} proteins, {len(gmst_scores):,} PSAURON scores")

    # --- 3-way merge ---
    all_genes = sorted(set(
        list(td2_seqs.keys()) +
        list(metaeuk_seqs.keys()) +
        list(gmst_seqs.keys())
    ))

    # Counters
    source_counts = {}
    n_filtered = 0
    n_rescued = 0
    n_zero_expr_dropped = 0
    n_total = 0

    # Per-predictor stats
    stats_by_source = {
        "td2":     {"lengths": [], "psaurons": []},
        "metaeuk": {"lengths": [], "psaurons": []},
        "gmst":    {"lengths": [], "psaurons": []},
    }

    with open(f"{species}.faa", "w") as faa, \
         open("merge_map.tsv", "w") as mapf:

        mapf.write("gene_id\tsource\tprotein_length\tpsauron_score\tcompleteness\tmean_tpm\t"
                   "td2_length\ttd2_psauron\ttd2_compl\t"
                   "metaeuk_length\tmetaeuk_psauron\tmetaeuk_compl\t"
                   "gmst_length\tgmst_psauron\tgmst_compl\n")

        for gene_id in all_genes:
            gene_tpm = gene_tpms.get(gene_id, 0.0)

            # Drop genes with zero expression (if expression data available)
            if has_expression and gene_tpm == 0.0:
                n_zero_expr_dropped += 1
                continue

            # Collect candidates: (source, seq, psauron, completeness)
            td2_seq = td2_seqs.get(gene_id)
            td2_ps = td2_scores.get(gene_id, 0.0)
            td2_c = td2_compl.get(gene_id, "unknown")
            td2_len = len(td2_seq) if td2_seq else 0

            met_seq = metaeuk_seqs.get(gene_id)
            met_ps = metaeuk_scores.get(gene_id, 0.0)
            met_c = classify_completeness(met_seq) if met_seq else "unknown"
            met_len = len(met_seq) if met_seq else 0

            gm_seq = gmst_seqs.get(gene_id)
            gm_ps = gmst_scores.get(gene_id, 0.0)
            gm_c = gmst_compl.get(gene_id, "unknown")
            gm_len = len(gm_seq) if gm_seq else 0

            # Standard filtering: PSAURON >= threshold
            candidates = []
            if td2_seq and td2_ps >= min_ps:
                candidates.append(("td2", td2_seq, td2_ps, td2_c))
            if met_seq and met_ps >= min_ps:
                candidates.append(("metaeuk", met_seq, met_ps, met_c))
            if gm_seq and gm_ps >= min_ps:
                candidates.append(("gmst", gm_seq, gm_ps, gm_c))

            rescued = False
            if not candidates:
                # Expression rescue: accept best candidate if gene has TPM >= threshold
                if has_expression and gene_tpm >= RESCUE_TPM_THRESHOLD:
                    all_cands = []
                    if td2_seq:
                        all_cands.append(("td2", td2_seq, td2_ps, td2_c))
                    if met_seq:
                        all_cands.append(("metaeuk", met_seq, met_ps, met_c))
                    if gm_seq:
                        all_cands.append(("gmst", gm_seq, gm_ps, gm_c))
                    if all_cands:
                        candidates = all_cands
                        rescued = True
                        n_rescued += 1

            if not candidates:
                n_filtered += 1
                continue

            # Rank: completeness → PSAURON → expression (mean TPM) → length
            candidates.sort(
                key=lambda x: (COMPL_RANK.get(x[3], 0), x[2], gene_tpm, len(x[1])),
                reverse=True
            )
            winner_src, winner_seq, winner_ps, winner_c = candidates[0]

            # Determine source label
            n_passing = len(candidates)
            if n_passing == 1:
                source_label = winner_src
            else:
                source_label = f"{'_'.join(sorted(c[0] for c in candidates))}_{winner_src}"

            if rescued:
                source_label += "_rescued"

            source_counts[source_label] = source_counts.get(source_label, 0) + 1

            faa.write(f">{gene_id}\n{winner_seq}\n")
            mapf.write(f"{gene_id}\t{source_label}\t{len(winner_seq)}\t{winner_ps:.4f}\t{winner_c}\t"
                       f"{gene_tpm:.4f}\t"
                       f"{td2_len}\t{td2_ps:.4f}\t{td2_c}\t"
                       f"{met_len}\t{met_ps:.4f}\t{met_c}\t"
                       f"{gm_len}\t{gm_ps:.4f}\t{gm_c}\n")

            # Accumulate per-predictor stats (for genes that had a prediction)
            if td2_seq:
                stats_by_source["td2"]["lengths"].append(td2_len)
                stats_by_source["td2"]["psaurons"].append(td2_ps)
            if met_seq:
                stats_by_source["metaeuk"]["lengths"].append(met_len)
                stats_by_source["metaeuk"]["psaurons"].append(met_ps)
            if gm_seq:
                stats_by_source["gmst"]["lengths"].append(gm_len)
                stats_by_source["gmst"]["psaurons"].append(gm_ps)

            n_total += 1

    # --- Statistics ---
    # Count genes where each predictor was the sole provider
    n_td2_only     = sum(v for k, v in source_counts.items() if k.replace("_rescued", "") == "td2")
    n_metaeuk_only = sum(v for k, v in source_counts.items() if k.replace("_rescued", "") == "metaeuk")
    n_gmst_only    = sum(v for k, v in source_counts.items() if k.replace("_rescued", "") == "gmst")
    n_multi = n_total - n_td2_only - n_metaeuk_only - n_gmst_only

    stats = f"""Merge statistics (min_psauron={min_ps})
Genes with TD2 prediction:       {len(td2_seqs):>10,}
Genes with MetaEuk prediction:   {len(metaeuk_seqs):>10,}
Genes with GeneMarkS-T prediction: {len(gmst_seqs):>10,}
Genes in union:                  {len(all_genes):>10,}
"""

    if has_expression:
        stats += f"""
Expression data:
  Samples loaded:                {len([d for d in args.salmon_quants.split(',') if d.strip()]):>10}
  Genes with expression data:    {len(gene_tpms):>10,}
  Zero-expression genes dropped: {n_zero_expr_dropped:>10,}
  Expression-rescued genes:      {n_rescued:>10,}
"""

    stats += f"""
After PSAURON filtering + completeness ranking:
  TD2-only genes:                {n_td2_only:>10,}
  MetaEuk-only genes:            {n_metaeuk_only:>10,}
  GeneMarkS-T-only genes:       {n_gmst_only:>10,}
  Multi-predictor genes:         {n_multi:>10,}
  Filtered (below threshold):   {n_filtered:>10,}
  Total merged proteins:         {n_total:>10,}

Per-source breakdown:
"""
    for label in sorted(source_counts, key=source_counts.get, reverse=True):
        cnt = source_counts[label]
        pct = cnt / n_total * 100 if n_total > 0 else 0.0
        stats += f"  {label:<40} {cnt:>8,}  ({pct:.1f}%)\n"

    stats += f"""
Per-predictor protein lengths:
  TD2:       mean={mean(stats_by_source['td2']['lengths']):.0f}  median={median(stats_by_source['td2']['lengths'])}
  MetaEuk:   mean={mean(stats_by_source['metaeuk']['lengths']):.0f}  median={median(stats_by_source['metaeuk']['lengths'])}
  GeneMarkST: mean={mean(stats_by_source['gmst']['lengths']):.0f}  median={median(stats_by_source['gmst']['lengths'])}

Per-predictor PSAURON scores:
  TD2:       mean={mean(stats_by_source['td2']['psaurons']):.3f}  median={median(stats_by_source['td2']['psaurons']):.3f}
  MetaEuk:   mean={mean(stats_by_source['metaeuk']['psaurons']):.3f}  median={median(stats_by_source['metaeuk']['psaurons']):.3f}
  GeneMarkST: mean={mean(stats_by_source['gmst']['psaurons']):.3f}  median={median(stats_by_source['gmst']['psaurons']):.3f}
"""

    with open("merge_stats.txt", "w") as f:
        f.write(stats)

    print(stats)


if __name__ == "__main__":
    main()
