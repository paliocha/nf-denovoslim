#!/usr/bin/env python3
"""Generate a thinning/slimming summary report.

Usage:
    thinning_report.py <species_label> \\
        <trinity_fasta> <supertranscripts_fasta> \\
        <grouper_clust> <orf_to_gene_map> <faa> \\
        <initial_quant_dirs> <final_quant_dirs> \\
        <busco_trinity_summary> <busco_final_summary> \\
        <id_validation> <sortmerna_logs> <taxonomy_breakdown>

Outputs:
    {species_label}_thinning_report.txt
"""

import csv
import os
import re
import sys
import statistics
from collections import defaultdict


def count_fasta(path):
    """Count sequences and collect lengths from a FASTA file."""
    lengths = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                lengths.append(0)
            elif lengths:
                lengths[-1] += len(line.strip())
    return lengths


def parse_sortmerna_log(path):
    """Parse a SortMeRNA v4 log file.

    Returns dict with keys: total, rrna, non_rrna, pct_rrna (or None on failure).
    SortMeRNA v4 logs contain lines like:
        Total reads = 12345678
        Total reads passing E-value threshold = 678901 (5.50%)
        Total reads failing E-value threshold = 11666777 (94.50%)
    """
    result = {}
    with open(path) as f:
        for line in f:
            m = re.search(r'Total reads = (\d+)', line)
            if m and 'passing' not in line and 'failing' not in line:
                result['total'] = int(m.group(1))
            m = re.search(r'Total reads passing E-value threshold = (\d+)', line)
            if m:
                result['rrna'] = int(m.group(1))
            m = re.search(r'Total reads failing E-value threshold = (\d+)', line)
            if m:
                result['non_rrna'] = int(m.group(1))
    if 'total' in result and 'rrna' in result:
        result['pct_rrna'] = result['rrna'] / result['total'] * 100 if result['total'] > 0 else 0.0
    return result if 'total' in result else None


def parse_quant_sf(path):
    """Parse a Salmon quant.sf file, return dict of Name -> (TPM, NumReads)."""
    quants = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            length = int(parts[1])
            tpm = float(parts[3])
            num_reads = float(parts[4])
            quants[name] = (tpm, num_reads)
    return quants


def aggregate_quants(quant_dirs, sf_name="quant.sf"):
    """Aggregate TPM and NumReads across samples, return per-gene means."""
    gene_tpm = defaultdict(list)
    gene_reads = defaultdict(list)
    n_samples = 0
    for qdir in quant_dirs:
        sf = os.path.join(qdir, sf_name)
        if not os.path.isfile(sf):
            continue
        n_samples += 1
        quants = parse_quant_sf(sf)
        for name, (tpm, nreads) in quants.items():
            gene_tpm[name].append(tpm)
            gene_reads[name].append(nreads)
    return gene_tpm, gene_reads, n_samples


def fmt(number):
    """Format a number with comma separators."""
    if isinstance(number, float):
        return f"{number:,.2f}"
    return f"{number:,}"


def main():
    if len(sys.argv) < 13:
        print(f"Usage: {sys.argv[0]} <species_label> <trinity.fasta> "
              f"<supertranscripts.fasta> <grouper_clust> <orf_map> <faa> "
              f"<initial_quant_dirs> <final_quant_dirs> <busco_trinity_summary> "
              f"<busco_final_summary> <id_validation> [sortmerna_logs] "
              f"[taxonomy_breakdown]",
              file=sys.stderr)
        sys.exit(1)

    species       = sys.argv[1]
    trinity_fasta = sys.argv[2]
    st_fasta      = sys.argv[3]
    grouper_clust = sys.argv[4]
    orf_map       = sys.argv[5]
    faa_file      = sys.argv[6]
    initial_qdir  = sys.argv[7]   # comma-separated list of dirs
    final_qdir    = sys.argv[8]   # comma-separated list of dirs
    busco_trinity = sys.argv[9]
    busco_final   = sys.argv[10]
    id_validation = sys.argv[11]
    sortmerna_arg = sys.argv[12] if len(sys.argv) > 12 else ""
    taxonomy_arg  = sys.argv[13] if len(sys.argv) > 13 else ""

    # ── SortMeRNA rRNA removal stats ───────────────────────────────────

    sortmerna_stats = []
    if sortmerna_arg:
        for logpath in sortmerna_arg.split(","):
            logpath = logpath.strip()
            if logpath and os.path.isfile(logpath):
                parsed = parse_sortmerna_log(logpath)
                if parsed:
                    parsed['sample'] = os.path.basename(os.path.dirname(logpath)) or logpath
                    sortmerna_stats.append(parsed)

    # ── Taxonomy breakdown ───────────────────────────────────────────────

    taxonomy_rows = []
    if taxonomy_arg and os.path.isfile(taxonomy_arg):
        with open(taxonomy_arg) as f:
            header = f.readline()  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    taxonomy_rows.append({
                        'category': parts[0],
                        'count': int(parts[1]),
                        'percent': parts[2],
                        'status': parts[3],
                    })

    # ── Sequence counts and lengths ──────────────────────────────────────

    trinity_lengths  = count_fasta(trinity_fasta)
    st_lengths       = count_fasta(st_fasta)
    protein_lengths  = count_fasta(faa_file)

    n_trinity  = len(trinity_lengths)
    n_genes    = len(st_lengths)
    n_proteins = len(protein_lengths)

    # ── Grouper cluster stats ────────────────────────────────────────────

    cluster_sizes = defaultdict(int)
    with open(grouper_clust) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                cluster_sizes[parts[1]] += 1
    cluster_counts = list(cluster_sizes.values())

    # ── ORF map stats ────────────────────────────────────────────────────

    psauron_scores = []
    orf_lengths = []
    with open(orf_map) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                psauron_scores.append(float(parts[2]))
                orf_lengths.append(int(parts[3]))

    # ── Expression stats (initial = transcript-level) ────────────────────

    initial_dirs = [d for d in initial_qdir.split(",") if d.strip()]
    init_tpm, init_reads, n_samples_init = aggregate_quants(initial_dirs)

    # Per-transcript mean TPM across samples
    init_mean_tpms = [statistics.mean(v) for v in init_tpm.values() if v]

    # ── Expression stats (final = gene-level SuperTranscripts) ───────────

    final_dirs = [d for d in final_qdir.split(",") if d.strip()]
    final_tpm, final_reads, n_samples_final = aggregate_quants(final_dirs)

    final_mean_tpms = [statistics.mean(v) for v in final_tpm.values() if v]

    # Zero-expression genes (mean TPM == 0 across all samples)
    n_zero_init  = sum(1 for t in init_mean_tpms if t == 0)
    n_zero_final = sum(1 for t in final_mean_tpms if t == 0)

    # ── BUSCO summaries ────────────────────────────────────────────────

    busco_trinity_text = ""
    if os.path.isfile(busco_trinity):
        with open(busco_trinity) as f:
            busco_trinity_text = f.read().strip()

    busco_final_text = ""
    if os.path.isfile(busco_final):
        with open(busco_final) as f:
            busco_final_text = f.read().strip()

    # ── ID validation ────────────────────────────────────────────────────

    id_val_text = ""
    if os.path.isfile(id_validation):
        with open(id_validation) as f:
            id_val_text = f.read().strip()

    # ── Write report ─────────────────────────────────────────────────────

    out = f"{species}_thinning_report.txt"
    with open(out, "w") as r:
        r.write(f"{'=' * 72}\n")
        r.write(f"  nf-denovoslim Thinning Report — {species}\n")
        r.write(f"{'=' * 72}\n\n")

        # --- Section 1: rRNA removal ---
        r.write("1. rRNA REMOVAL (SortMeRNA)\n")
        r.write("-" * 40 + "\n")
        if sortmerna_stats:
            total_reads = sum(s['total'] for s in sortmerna_stats)
            total_rrna  = sum(s['rrna'] for s in sortmerna_stats)
            total_kept  = sum(s.get('non_rrna', s['total'] - s['rrna']) for s in sortmerna_stats)
            overall_pct = total_rrna / total_reads * 100 if total_reads > 0 else 0.0

            r.write(f"  Samples processed:    {len(sortmerna_stats)}\n")
            r.write(f"  Total reads:          {fmt(total_reads)}\n")
            r.write(f"  rRNA reads:           {fmt(total_rrna)}  ({overall_pct:.2f}%)\n")
            r.write(f"  Non-rRNA (kept):      {fmt(total_kept)}  ({100 - overall_pct:.2f}%)\n\n")

            pct_values = [s['pct_rrna'] for s in sortmerna_stats]
            r.write(f"  Per-sample rRNA %:    min {min(pct_values):.2f}%,  "
                    f"max {max(pct_values):.2f}%,  "
                    f"mean {statistics.mean(pct_values):.2f}%\n")
        else:
            r.write("  (SortMeRNA log data not available)\n")
        r.write("\n")

        # --- Section 2: Assembly reduction ---
        r.write("2. ASSEMBLY REDUCTION\n")
        r.write("-" * 40 + "\n")
        r.write(f"  Original Trinity transcripts:     {fmt(n_trinity)}\n")
        r.write(f"  After Corset -> Lace:              {fmt(n_genes)}  "
                f"({n_genes/n_trinity*100:.1f}% of original)\n")
        r.write(f"  Genes with predicted ORF:          {fmt(n_proteins)}  "
                f"({n_proteins/n_genes*100:.1f}% of genes)\n")
        r.write(f"\n  Overall reduction: {n_trinity} -> {n_genes} "
                f"({n_trinity/n_genes:.1f}x collapse)\n\n")

        # --- Section 3: Taxonomy filter ---
        r.write("3. TAXONOMY FILTER (MMseqs2)\n")
        r.write("-" * 40 + "\n")
        if taxonomy_rows:
            # Find TOTAL row for summary
            total_row = next((row for row in taxonomy_rows if row['category'] == 'TOTAL'), None)
            kept_rows = [row for row in taxonomy_rows if row['status'] == 'KEPT']
            removed_rows = [row for row in taxonomy_rows if row['status'] == 'REMOVED']
            kept_count = sum(row['count'] for row in kept_rows)
            removed_count = sum(row['count'] for row in removed_rows)

            if total_row:
                total_count = total_row['count']
                r.write(f"  SuperTranscripts classified:  {fmt(total_count)}\n")
                r.write(f"  Kept (Viridiplantae + no-hit): {fmt(kept_count)}  "
                        f"({kept_count/total_count*100:.1f}%)\n")
                r.write(f"  Removed (non-plant):           {fmt(removed_count)}  "
                        f"({removed_count/total_count*100:.1f}%)\n")
            r.write("\n")

            # Breakdown table
            cat_width = max(len(row['category']) for row in taxonomy_rows) + 2
            r.write(f"  {'Category':<{cat_width}} {'Count':>10} {'Percent':>10}  Status\n")
            r.write(f"  {'-' * cat_width} {'-' * 10} {'-' * 10}  {'-' * 7}\n")
            for row in taxonomy_rows:
                if row['category'] == 'TOTAL':
                    r.write(f"  {'-' * cat_width} {'-' * 10} {'-' * 10}  {'-' * 7}\n")
                r.write(f"  {row['category']:<{cat_width}} {row['count']:>10,} "
                        f"{row['percent']:>10}  {row['status']}\n")
        else:
            r.write("  (taxonomy breakdown not available)\n")
        r.write("\n")

        # --- Section 4: Sequence lengths ---
        r.write("4. SEQUENCE LENGTH DISTRIBUTION (nt)\n")
        r.write("-" * 40 + "\n")
        r.write(f"  {'Stage':<30} {'N':>8} {'Mean':>10} {'Median':>10} "
                f"{'Min':>8} {'Max':>10}\n")
        for label, lens in [("Trinity transcripts", trinity_lengths),
                            ("SuperTranscripts", st_lengths)]:
            if lens:
                r.write(f"  {label:<30} {len(lens):>8,} {statistics.mean(lens):>10,.1f} "
                        f"{statistics.median(lens):>10,.1f} "
                        f"{min(lens):>8,} {max(lens):>10,}\n")

        # Length increase from Trinity to SuperTranscripts
        if trinity_lengths and st_lengths:
            trinity_median = statistics.median(trinity_lengths)
            st_median = statistics.median(st_lengths)
            median_increase = (st_median - trinity_median) / trinity_median * 100
            r.write(f"\n  Median length increase (Trinity → LACE): "
                    f"{trinity_median:,.0f} bp → {st_median:,.0f} bp "
                    f"(+{median_increase:.1f}%, {st_median/trinity_median:.1f}x)\n")
        r.write("\n")

        # --- Section 5: Protein/ORF stats ---
        r.write("5. PREDICTED PROTEINS\n")
        r.write("-" * 40 + "\n")
        r.write(f"  Total proteins:      {fmt(n_proteins)}\n")
        if orf_lengths:
            r.write(f"  Mean protein length: {statistics.mean(orf_lengths):,.1f} aa\n")
            r.write(f"  Median protein len:  {statistics.median(orf_lengths):,.1f} aa\n")
        if psauron_scores:
            r.write(f"  Mean PSAURON score:  {statistics.mean(psauron_scores):.3f}\n")
            r.write(f"  Median PSAURON:      {statistics.median(psauron_scores):.3f}\n")
        r.write("\n")

        # --- Section 6: Corset cluster sizes ---
        r.write("6. CORSET CLUSTER SIZES\n")
        r.write("-" * 40 + "\n")
        if cluster_counts:
            r.write(f"  Total clusters (genes): {fmt(len(cluster_counts))}\n")
            r.write(f"  Mean transcripts/gene:  {statistics.mean(cluster_counts):.2f}\n")
            r.write(f"  Median transcripts/gene: {statistics.median(cluster_counts):.1f}\n")
            r.write(f"  Max transcripts/gene:   {max(cluster_counts)}\n")
            singleton = sum(1 for c in cluster_counts if c == 1)
            r.write(f"  Singleton genes:        {fmt(singleton)} "
                    f"({singleton/len(cluster_counts)*100:.1f}%)\n")
        r.write("\n")

        # --- Section 7: Expression summary ---
        r.write("7. EXPRESSION SUMMARY (mean TPM across samples)\n")
        r.write("-" * 40 + "\n")
        r.write(f"  {'Metric':<35} {'Transcript':>12} {'Gene':>12}\n")
        if init_mean_tpms and final_mean_tpms:
            r.write(f"  {'Samples quantified':<35} {n_samples_init:>12} "
                    f"{n_samples_final:>12}\n")
            r.write(f"  {'Features quantified':<35} {len(init_mean_tpms):>12,} "
                    f"{len(final_mean_tpms):>12,}\n")
            r.write(f"  {'Mean TPM (across features)':<35} "
                    f"{statistics.mean(init_mean_tpms):>12.2f} "
                    f"{statistics.mean(final_mean_tpms):>12.2f}\n")
            r.write(f"  {'Median TPM':<35} "
                    f"{statistics.median(init_mean_tpms):>12.2f} "
                    f"{statistics.median(final_mean_tpms):>12.2f}\n")
            nz_init = [t for t in init_mean_tpms if t > 0]
            nz_final = [t for t in final_mean_tpms if t > 0]
            if nz_init and nz_final:
                r.write(f"  {'Median TPM (non-zero only)':<35} "
                        f"{statistics.median(nz_init):>12.2f} "
                        f"{statistics.median(nz_final):>12.2f}\n")
            r.write(f"  {'Zero-expression features':<35} "
                    f"{n_zero_init:>12,} {n_zero_final:>12,}\n")
        else:
            r.write("  (expression data not available)\n")
        r.write("\n")

        # --- Section 8: BUSCO ---
        r.write("8. BUSCO COMPLETENESS\n")
        r.write("-" * 40 + "\n")

        r.write("  --- Trinity baseline (transcriptome mode) ---\n")
        if busco_trinity_text:
            for line in busco_trinity_text.splitlines():
                r.write(f"  {line}\n")
        else:
            r.write("  (BUSCO Trinity summary not available)\n")
        r.write("\n")

        r.write("  --- Final proteins (protein mode) ---\n")
        if busco_final_text:
            for line in busco_final_text.splitlines():
                r.write(f"  {line}\n")
        else:
            r.write("  (BUSCO final summary not available)\n")
        r.write("\n")

        # --- Section 9: ID validation ---
        r.write("9. ID CONSISTENCY CHECK\n")
        r.write("-" * 40 + "\n")
        if id_val_text:
            for line in id_val_text.splitlines():
                r.write(f"  {line}\n")
        else:
            r.write("  (validation data not available)\n")
        r.write("\n")

        r.write(f"{'=' * 72}\n")
        r.write(f"  Report generated by nf-denovoslim\n")
        r.write(f"{'=' * 72}\n")

    print(f"Report written to {out}")


if __name__ == "__main__":
    main()
