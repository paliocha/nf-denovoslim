#!/usr/bin/env python3
"""Generate a thinning/slimming summary report.

Usage:
    thinning_report.py --species BMAX --trinity trinity.fasta \\
        --representatives reps.fasta --clusters corset.clust \\
        --merge-map map.tsv --faa proteins.faa \\
        --initial-quants dir1,dir2 --final-quants dir1,dir2 \\
        --busco-trinity short_summary.txt --busco-final short_summary.txt \\
        --id-validation report.txt [--sortmerna-logs log1,log2] \\
        [--taxonomy breakdown.tsv] [--transannot annot.tsv] \\
        [--merge-stats stats.txt] [--dedup-stats stats.txt] \\
        [--protein-dedup-stats stats.txt]

Outputs:
    {species}_thinning_report.txt
"""

import argparse
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

    Handles two formats:
      1. Summary format (non-verbose):
           Total reads = 12345678
           Total reads passing E-value threshold = 678901 (5.50%)
           Total reads failing E-value threshold = 11666777 (94.50%)
      2. Verbose (-v) format:
           [count_reads:945] done count. ... Total reads: 64536010
           [align2:133] Processor N ... Aligned reads (passing E-value): M
    """
    result = {}
    v_total = None
    v_aligned_sum = 0
    with open(path) as f:
        for line in f:
            # --- Summary format (non-verbose) ---
            m = re.search(r'Total reads = (\d+)', line)
            if m and 'passing' not in line and 'failing' not in line:
                result['total'] = int(m.group(1))
            m = re.search(r'Total reads passing E-value threshold = (\d+)', line)
            if m:
                result['rrna'] = int(m.group(1))
            m = re.search(r'Total reads failing E-value threshold = (\d+)', line)
            if m:
                result['non_rrna'] = int(m.group(1))
            # --- Verbose format ---
            m = re.search(r'\[count_reads:\d+\] done count\..+Total reads: (\d+)', line)
            if m:
                v_total = int(m.group(1))
            m = re.search(
                r'\[align2:\d+\] Processor \d+.*Aligned reads \(passing E-value\): (\d+)',
                line,
            )
            if m:
                v_aligned_sum += int(m.group(1))
    # Prefer summary format if available
    if 'total' in result and 'rrna' in result:
        result['pct_rrna'] = result['rrna'] / result['total'] * 100 if result['total'] > 0 else 0.0
        return result
    # Fall back to verbose format
    if v_total is not None and v_total > 0:
        result['total'] = v_total
        result['rrna'] = v_aligned_sum
        result['non_rrna'] = v_total - v_aligned_sum
        result['pct_rrna'] = v_aligned_sum / v_total * 100
        return result
    return None


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
    parser = argparse.ArgumentParser(
        description='Generate nf-denovoslim thinning report')
    parser.add_argument('--species', required=True,
                        help='Species label (e.g. BMAX)')
    parser.add_argument('--trinity', required=True,
                        help='Trinity FASTA')
    parser.add_argument('--representatives', required=True,
                        help='Representative sequences FASTA')
    parser.add_argument('--clusters', required=True,
                        help='Corset cluster file')
    parser.add_argument('--merge-map', required=True,
                        help='Merge prediction map TSV')
    parser.add_argument('--faa', required=True,
                        help='Final protein FASTA')
    parser.add_argument('--initial-quants', required=True,
                        help='Comma-separated initial quant dirs')
    parser.add_argument('--final-quants', required=True,
                        help='Comma-separated final quant dirs')
    parser.add_argument('--busco-trinity', required=True,
                        help='BUSCO Trinity summary file')
    parser.add_argument('--busco-reps', default='',
                        help='BUSCO representatives summary file')
    parser.add_argument('--busco-corrected', default='',
                        help='BUSCO corrected representatives summary file')
    parser.add_argument('--busco-final', required=True,
                        help='BUSCO final summary file')
    parser.add_argument('--id-validation', required=True,
                        help='ID validation report')
    parser.add_argument('--sortmerna-logs', default='',
                        help='Comma-separated SortMeRNA log files')
    parser.add_argument('--taxonomy', default='',
                        help='Taxonomy breakdown TSV')
    parser.add_argument('--transannot', default='',
                        help='TransAnnot annotation TSV')
    parser.add_argument('--merge-stats', default='',
                        help='Merge statistics text file')
    parser.add_argument('--dedup-stats', default='',
                        help='Nucleotide dedup statistics text file')
    parser.add_argument('--protein-dedup-stats', default='',
                        help='Protein dedup statistics text file')
    args = parser.parse_args()

    species       = args.species
    trinity_fasta = args.trinity
    rep_fasta     = args.representatives
    grouper_clust = args.clusters
    merge_map     = args.merge_map
    faa_file      = args.faa
    initial_qdir  = args.initial_quants
    final_qdir    = args.final_quants
    busco_trinity = args.busco_trinity
    busco_reps    = args.busco_reps
    busco_corrected = args.busco_corrected
    busco_final   = args.busco_final
    id_validation = args.id_validation
    sortmerna_arg = args.sortmerna_logs
    taxonomy_arg  = args.taxonomy
    annot_arg     = args.transannot
    merge_stats   = args.merge_stats
    dedup_stats   = args.dedup_stats
    protein_dedup_stats_arg = args.protein_dedup_stats

    # -- SortMeRNA rRNA removal stats --

    sortmerna_stats = []
    if sortmerna_arg:
        for logpath in sortmerna_arg.split(","):
            logpath = logpath.strip()
            if logpath and os.path.isfile(logpath):
                parsed = parse_sortmerna_log(logpath)
                if parsed:
                    parsed['sample'] = os.path.basename(os.path.dirname(logpath)) or logpath
                    sortmerna_stats.append(parsed)

    # -- Taxonomy breakdown --

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

    # -- TransAnnot functional annotation --

    annot_by_protein = defaultdict(set)
    pfam_counts = defaultdict(int)
    if annot_arg and os.path.isfile(annot_arg):
        with open(annot_arg) as f:
            for line in f:
                if line.startswith("queryID"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 10:
                    query_id = parts[0]
                    target_desc = parts[4]
                    db_name = parts[9].strip()
                    annot_by_protein[query_id].add(db_name)
                    if 'pfam' in db_name.lower():
                        pfam_fam = target_desc.split()[0] if target_desc else 'unknown'
                        pfam_counts[pfam_fam] += 1

    # -- Merge stats (pre-formatted text) --

    merge_stats_text = ""
    if merge_stats and os.path.isfile(merge_stats):
        with open(merge_stats) as f:
            merge_stats_text = f.read().strip()

    # -- Dedup stats (pre-formatted text) --

    dedup_stats_text = ""
    if dedup_stats and os.path.isfile(dedup_stats):
        with open(dedup_stats) as f:
            dedup_stats_text = f.read().strip()

    # -- Protein dedup stats (pre-formatted text) --

    protein_dedup_stats_text = ""
    if protein_dedup_stats_arg and os.path.isfile(protein_dedup_stats_arg):
        with open(protein_dedup_stats_arg) as f:
            protein_dedup_stats_text = f.read().strip()

    # -- Sequence counts and lengths --

    trinity_lengths  = count_fasta(trinity_fasta)
    rep_lengths      = count_fasta(rep_fasta)
    protein_lengths  = count_fasta(faa_file)

    n_trinity  = len(trinity_lengths)
    n_reps     = len(rep_lengths)
    n_proteins = len(protein_lengths)

    # -- Grouper cluster stats --

    cluster_sizes = defaultdict(int)
    with open(grouper_clust) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                cluster_sizes[parts[1]] += 1
    cluster_counts = list(cluster_sizes.values())

    # -- Merge map stats (TD2 + MetaEuk breakdown) --

    merge_sources = defaultdict(int)
    merge_psaurons = []
    merge_lengths = []
    merge_completeness = defaultdict(int)
    td2_psaurons = []
    metaeuk_psaurons = []
    gmst_psaurons = []
    td2_lengths = []
    metaeuk_lengths = []
    gmst_lengths = []

    merge_mean_tpms = []   # per-gene mean TPM from merge step (if available)

    if os.path.isfile(merge_map):
        with open(merge_map) as f:
            header_line = f.readline().strip()
            cols = header_line.split("\t")
            # Build column index map for robust parsing
            ci = {name: i for i, name in enumerate(cols)}
            has_tpm = "mean_tpm" in ci
            has_3way = "gmst_length" in ci

            for line in f:
                parts = line.strip().split("\t")
                if has_3way and len(parts) >= len(cols):
                    source = parts[ci["source"]]
                    prot_len = int(parts[ci["protein_length"]])
                    psauron = float(parts[ci["psauron_score"]])
                    compl = parts[ci["completeness"]]
                    td2_len = int(parts[ci["td2_length"]])
                    td2_ps = float(parts[ci["td2_psauron"]])
                    met_len = int(parts[ci["metaeuk_length"]])
                    met_ps = float(parts[ci["metaeuk_psauron"]])
                    gm_len = int(parts[ci["gmst_length"]])
                    gm_ps = float(parts[ci["gmst_psauron"]])

                    merge_sources[source] += 1
                    merge_psaurons.append(psauron)
                    merge_lengths.append(prot_len)
                    merge_completeness[compl] += 1

                    if has_tpm:
                        merge_mean_tpms.append(float(parts[ci["mean_tpm"]]))

                    if td2_len > 0:
                        td2_psaurons.append(td2_ps)
                        td2_lengths.append(td2_len)
                    if met_len > 0:
                        metaeuk_psaurons.append(met_ps)
                        metaeuk_lengths.append(met_len)
                    if gm_len > 0:
                        gmst_psaurons.append(gm_ps)
                        gmst_lengths.append(gm_len)
                elif len(parts) >= 8:  # Legacy 2-way format
                    source = parts[1]
                    prot_len = int(parts[2])
                    psauron = float(parts[3])
                    td2_len = int(parts[4])
                    td2_ps = float(parts[5])
                    met_len = int(parts[6])
                    met_ps = float(parts[7])

                    merge_sources[source] += 1
                    merge_psaurons.append(psauron)
                    merge_lengths.append(prot_len)

                    if td2_len > 0:
                        td2_psaurons.append(td2_ps)
                        td2_lengths.append(td2_len)
                    if met_len > 0:
                        metaeuk_psaurons.append(met_ps)
                        metaeuk_lengths.append(met_len)

    # -- Expression stats (initial = transcript-level) --

    initial_dirs = [d for d in initial_qdir.split(",") if d.strip()]
    init_tpm, init_reads, n_samples_init = aggregate_quants(initial_dirs)
    init_mean_tpms = [statistics.mean(v) for v in init_tpm.values() if v]

    # -- Expression stats (final = gene-level representatives) --

    final_dirs = [d for d in final_qdir.split(",") if d.strip()]
    final_tpm, final_reads, n_samples_final = aggregate_quants(final_dirs)
    final_mean_tpms = [statistics.mean(v) for v in final_tpm.values() if v]

    n_zero_init  = sum(1 for t in init_mean_tpms if t == 0)
    n_zero_final = sum(1 for t in final_mean_tpms if t == 0)

    # -- BUSCO summaries --

    busco_trinity_text = ""
    if os.path.isfile(busco_trinity):
        with open(busco_trinity) as f:
            busco_trinity_text = f.read().strip()

    busco_reps_text = ""
    if busco_reps and os.path.isfile(busco_reps):
        with open(busco_reps) as f:
            busco_reps_text = f.read().strip()

    busco_corrected_text = ""
    if busco_corrected and os.path.isfile(busco_corrected):
        with open(busco_corrected) as f:
            busco_corrected_text = f.read().strip()

    busco_final_text = ""
    if os.path.isfile(busco_final):
        with open(busco_final) as f:
            busco_final_text = f.read().strip()

    # -- ID validation --

    id_val_text = ""
    if os.path.isfile(id_validation):
        with open(id_validation) as f:
            id_val_text = f.read().strip()

    # -- Write report --

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
        r.write(f"  Corset clusters (genes):           {fmt(len(cluster_counts))}\n")
        # Parse taxonomy-filter input count from dedup stats
        if dedup_stats_text and 'Input sequences:' in dedup_stats_text:
            n_after_tax = dedup_stats_text.split('Input sequences:')[1].split()[0].strip()
            r.write(f"  After taxonomy filter:             {n_after_tax}\n")
        r.write(f"  After nucleotide dedup:            {fmt(n_reps)}\n")
        if dedup_stats_text:
            for line in dedup_stats_text.strip().splitlines():
                r.write(f"    {line.strip()}\n")
        r.write(f"  Genes with predicted protein:      {fmt(n_proteins)}  "
                f"({n_proteins/n_reps*100:.1f}% of representatives)\n")
        r.write(f"\n  Overall reduction: {n_trinity} -> {n_reps} "
                f"({n_trinity/n_reps:.1f}x collapse)\n\n")

        # --- Section 3: Taxonomy filter ---
        r.write("3. TAXONOMY FILTER (MMseqs2)\n")
        r.write("-" * 40 + "\n")
        if taxonomy_rows:
            total_row = next((row for row in taxonomy_rows if row['category'] == 'TOTAL'), None)
            kept_rows = [row for row in taxonomy_rows if row['status'] == 'KEPT']
            removed_rows = [row for row in taxonomy_rows if row['status'] == 'REMOVED']
            kept_count = sum(row['count'] for row in kept_rows)
            removed_count = sum(row['count'] for row in removed_rows)

            if total_row:
                total_count = total_row['count']
                r.write(f"  Representatives classified:   {fmt(total_count)}\n")
                r.write(f"  Kept (Viridiplantae + no-hit): {fmt(kept_count)}  "
                        f"({kept_count/total_count*100:.1f}%)\n")
                r.write(f"  Removed (non-plant):           {fmt(removed_count)}  "
                        f"({removed_count/total_count*100:.1f}%)\n")
            r.write("\n")

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
                            ("Representatives", rep_lengths)]:
            if lens:
                r.write(f"  {label:<30} {len(lens):>8,} {statistics.mean(lens):>10,.1f} "
                        f"{statistics.median(lens):>10,.1f} "
                        f"{min(lens):>8,} {max(lens):>10,}\n")
        r.write("\n")

        # --- Section 5: Predicted proteins ---
        r.write("5. PREDICTED PROTEINS (TD2 + MetaEuk + GeneMarkS-T merge)\n")
        r.write("-" * 40 + "\n")
        r.write(f"  Total merged proteins: {fmt(sum(merge_sources.values()))}\n")
        r.write(f"  After protein dedup:   {fmt(n_proteins)}\n")
        if merge_lengths:
            r.write(f"  Mean protein length:   {statistics.mean(merge_lengths):,.1f} aa\n")
            r.write(f"  Median protein len:    {statistics.median(merge_lengths):,.1f} aa\n")
        if merge_psaurons:
            r.write(f"  Mean PSAURON score:    {statistics.mean(merge_psaurons):.3f}\n")
            r.write(f"  Median PSAURON:        {statistics.median(merge_psaurons):.3f}\n")
        if merge_mean_tpms:
            n_rescued = sum(1 for s in merge_sources if "rescued" in s)
            n_rescued_genes = sum(v for k, v in merge_sources.items() if "rescued" in k)
            r.write(f"  Expression-rescued genes: {fmt(n_rescued_genes)}\n")
            r.write(f"  Mean gene TPM (merged): {statistics.mean(merge_mean_tpms):.2f}\n")
            r.write(f"  Median gene TPM:        {statistics.median(merge_mean_tpms):.2f}\n")
        r.write("\n")

        # Completeness breakdown
        if merge_completeness:
            r.write("  Completeness breakdown (best-per-gene):\n")
            for compl_type in sorted(merge_completeness, key=merge_completeness.get, reverse=True):
                cnt = merge_completeness[compl_type]
                total_merged = sum(merge_sources.values())
                pct = cnt / total_merged * 100 if total_merged > 0 else 0.0
                r.write(f"    {compl_type:<25} {cnt:>10,}  ({pct:.1f}%)\n")
            r.write("\n")

        # Source breakdown (all labels from 3-way merge)
        if merge_sources:
            r.write("  Prediction source breakdown:\n")
            total_merged = sum(merge_sources.values())
            for source in sorted(merge_sources, key=merge_sources.get, reverse=True):
                cnt = merge_sources[source]
                pct = cnt / total_merged * 100 if total_merged > 0 else 0.0
                r.write(f"    {source:<40} {cnt:>10,}  ({pct:.1f}%)\n")
            r.write("\n")

        # Per-predictor stats (3-way)
        predictor_data = [
            ('TD2', td2_lengths, td2_psaurons),
            ('MetaEuk', metaeuk_lengths, metaeuk_psaurons),
            ('GeneMarkS-T', gmst_lengths, gmst_psaurons),
        ]
        active = [(name, lens, pss) for name, lens, pss in predictor_data if lens]
        if active:
            col_w = 14
            r.write("  Per-predictor statistics:\n")
            header_line = f"    {'Metric':<30}"
            for name, _, _ in active:
                header_line += f" {name:>{col_w}}"
            r.write(header_line + "\n")
            for metric, getter in [
                ('Genes with prediction', lambda L, _: f"{len(L):>,}"),
                ('Mean protein length (aa)', lambda L, _: f"{statistics.mean(L):>,.1f}" if L else 'N/A'),
                ('Median protein length', lambda L, _: f"{statistics.median(L):>,.1f}" if L else 'N/A'),
                ('Mean PSAURON score', lambda _, P: f"{statistics.mean(P):>.3f}" if P else 'N/A'),
                ('Median PSAURON', lambda _, P: f"{statistics.median(P):>.3f}" if P else 'N/A'),
            ]:
                row_line = f"    {metric:<30}"
                for _, lens, pss in active:
                    row_line += f" {getter(lens, pss):>{col_w}}"
                r.write(row_line + "\n")
        r.write("\n")

        # Protein dedup stats
        if protein_dedup_stats_text:
            r.write("  Protein-level dedup (95% aa identity):\n")
            for line in protein_dedup_stats_text.strip().splitlines():
                r.write(f"    {line.strip()}\n")
            r.write("\n")

        # Full merge stats text
        if merge_stats_text:
            r.write("  Full merge statistics:\n")
            for line in merge_stats_text.splitlines():
                r.write(f"    {line}\n")
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

        r.write("  --- Representatives (transcriptome mode) ---\n")
        if busco_reps_text:
            for line in busco_reps_text.splitlines():
                r.write(f"  {line}\n")
        else:
            r.write("  (BUSCO representatives summary not available)\n")
        r.write("\n")

        r.write("  --- Corrected representatives (transcriptome mode) ---\n")
        if busco_corrected_text:
            for line in busco_corrected_text.splitlines():
                r.write(f"  {line}\n")
        else:
            r.write("  (BUSCO corrected representatives summary not available)\n")
        r.write("\n")

        r.write("  --- Final proteins (protein mode, TD2+MetaEuk+GeneMarkS-T merged) ---\n")
        if busco_final_text:
            for line in busco_final_text.splitlines():
                r.write(f"  {line}\n")
        else:
            r.write("  (BUSCO final summary not available)\n")
        r.write("\n")

        # --- Section 9: Functional annotation ---
        r.write("9. FUNCTIONAL ANNOTATION (TransAnnot)\n")
        r.write("-" * 40 + "\n")
        if annot_by_protein:
            n_annotated = len(annot_by_protein)
            pct_annotated = n_annotated / n_proteins * 100 if n_proteins > 0 else 0.0

            all_dbs = set()
            for dbs in annot_by_protein.values():
                all_dbs.update(dbs)

            db_protein_count = defaultdict(int)
            for dbs in annot_by_protein.values():
                for db in dbs:
                    db_protein_count[db] += 1

            r.write(f"  Proteins with any annotation: {fmt(n_annotated)} / "
                    f"{fmt(n_proteins)}  ({pct_annotated:.1f}%)\n\n")

            r.write(f"  {'Database':<30} {'Proteins':>10} {'% of total':>12}\n")
            r.write(f"  {'-' * 30} {'-' * 10} {'-' * 12}\n")
            for db in sorted(db_protein_count, key=db_protein_count.get, reverse=True):
                cnt = db_protein_count[db]
                pct = cnt / n_proteins * 100 if n_proteins > 0 else 0.0
                r.write(f"  {db:<30} {cnt:>10,} {pct:>11.1f}%\n")

            if pfam_counts:
                r.write(f"\n  Top 15 Pfam domains:\n")
                for domain, count in sorted(pfam_counts.items(),
                                            key=lambda x: x[1],
                                            reverse=True)[:15]:
                    r.write(f"    {domain:<40} {count:>8,}\n")
        else:
            r.write("  (annotation data not available)\n")
        r.write("\n")

        # --- Section 10: ID validation ---
        r.write("10. ID CONSISTENCY CHECK\n")
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
