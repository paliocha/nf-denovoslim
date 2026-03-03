#!/usr/bin/env python3
"""
filter_unmapped.py — Remove proteins from unmapped transcripts lacking expression.

For genome-guided runs, LOCUS_CLUSTER marks each transcript as gene-assigned,
intergenic, or unmapped.  Transcripts assigned to a gene or intergenic locus
already have evidence (genomic location).  Unmapped transcripts that DO have a
predicted ORF should still be retained if they show expression — but orphan
proteins from unmapped transcripts with no expression support are likely
noise.

Filtering rule for UNMAPPED transcripts:
  Keep if the transcript has TPM > min_tpm in >= min_samples samples.
  Remove otherwise.

All gene/intergenic transcripts are kept unconditionally.
If no locus map is provided, all proteins are kept (no-op).

Inputs:
  --proteins    FASTA of merged/extended proteins (one per gene)
  --locus-map   TSV from LOCUS_CLUSTER (transcript_id, gene_id, status)
  --quant-dirs  Salmon quant directories (one per sample, each with quant.sf)
  --min-tpm     Minimum TPM threshold [default: 1.0]
  --min-samples Minimum number of samples meeting TPM threshold [default: 2]
  --out         Output filtered protein FASTA
  --stats       Output stats file
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from Bio import SeqIO


def load_locus_map(map_file):
    """Load locus map, return set of unmapped transcript IDs."""
    unmapped = set()
    with open(map_file) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row['status'] == 'retained' and row.get('gene_id', '') == 'unmapped':
                unmapped.add(row['transcript_id'])
    return unmapped


def load_expression(quant_dirs, transcript_ids, min_tpm):
    """Load Salmon quant.sf files, count how many samples each transcript
    exceeds min_tpm.

    Returns dict of transcript_id -> number of samples with TPM >= min_tpm.
    """
    counts = defaultdict(int)
    for qdir in quant_dirs:
        sf_path = os.path.join(qdir, 'quant.sf')
        if not os.path.exists(sf_path):
            print(f"Warning: {sf_path} not found, skipping", file=sys.stderr)
            continue
        with open(sf_path) as fh:
            header = fh.readline()  # skip header
            for line in fh:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 4:
                    continue
                name = fields[0]
                tpm = float(fields[3])
                if name in transcript_ids and tpm >= min_tpm:
                    counts[name] += 1
    return counts


def protein_to_transcript(protein_id):
    """Extract transcript ID from protein header.

    Protein IDs follow patterns like:
      TRINITY_DN123_c0_s0_g1_i1.p1  ->  TRINITY_DN123_c0_s0_g1_i1
      gene_TRINITY_DN123_c0_s0_g1_i1  ->  TRINITY_DN123_c0_s0_g1_i1
    The transcript ID is the part before .p<N> suffix.
    """
    # Strip .p<N> suffix if present
    pid = protein_id
    if '.p' in pid:
        pid = pid.rsplit('.p', 1)[0]
    # Strip gene_ prefix if present (from metaeuk/gmst)
    if pid.startswith('gene_'):
        pid = pid[5:]
    return pid


def main():
    parser = argparse.ArgumentParser(
        description='Filter proteins from unmapped transcripts lacking expression')
    parser.add_argument('--proteins', required=True,
                        help='Input protein FASTA')
    parser.add_argument('--locus-map', required=True,
                        help='Locus map TSV from LOCUS_CLUSTER')
    parser.add_argument('--quant-dirs', nargs='+', required=True,
                        help='Salmon quant directories')
    parser.add_argument('--min-tpm', type=float, default=1.0,
                        help='Min TPM threshold [default: 1.0]')
    parser.add_argument('--min-samples', type=int, default=2,
                        help='Min samples exceeding TPM [default: 2]')
    parser.add_argument('--out', required=True,
                        help='Output filtered protein FASTA')
    parser.add_argument('--stats', required=True,
                        help='Output stats file')
    args = parser.parse_args()

    # Load locus map to identify unmapped transcripts
    unmapped_tx = load_locus_map(args.locus_map)
    print(f"Unmapped transcripts in locus map: {len(unmapped_tx)}",
          file=sys.stderr)

    # Read all protein IDs and map to transcript IDs
    proteins = list(SeqIO.parse(args.proteins, 'fasta'))
    print(f"Input proteins: {len(proteins)}", file=sys.stderr)

    # Find which proteins come from unmapped transcripts
    unmapped_proteins = {}  # protein_id -> transcript_id
    for rec in proteins:
        tx_id = protein_to_transcript(rec.id)
        if tx_id in unmapped_tx:
            unmapped_proteins[rec.id] = tx_id

    print(f"Proteins from unmapped transcripts: {len(unmapped_proteins)}",
          file=sys.stderr)

    if not unmapped_proteins:
        # No unmapped proteins — write everything through
        with open(args.out, 'w') as out_fh:
            SeqIO.write(proteins, out_fh, 'fasta')
        with open(args.stats, 'w') as sfh:
            sfh.write("Filter unmapped: no unmapped proteins found, all retained\n")
        print("No unmapped proteins, all retained", file=sys.stderr)
        return

    # Load expression data for unmapped transcripts only
    unmapped_tx_set = set(unmapped_proteins.values())
    expr_counts = load_expression(args.quant_dirs, unmapped_tx_set, args.min_tpm)

    # Decide which unmapped proteins to keep
    keep_proteins = set()
    remove_proteins = set()
    for pid, tx_id in unmapped_proteins.items():
        n_samples = expr_counts.get(tx_id, 0)
        if n_samples >= args.min_samples:
            keep_proteins.add(pid)
        else:
            remove_proteins.add(pid)

    print(f"Unmapped proteins kept (expressed): {len(keep_proteins)}",
          file=sys.stderr)
    print(f"Unmapped proteins removed: {len(remove_proteins)}",
          file=sys.stderr)

    # Write filtered FASTA
    n_written = 0
    with open(args.out, 'w') as out_fh:
        for rec in proteins:
            if rec.id not in remove_proteins:
                SeqIO.write(rec, out_fh, 'fasta')
                n_written += 1

    # Stats
    n_input = len(proteins)
    n_gene_intergenic = n_input - len(unmapped_proteins)
    n_removed = len(remove_proteins)
    n_kept_unmapped = len(keep_proteins)

    stats = (
        f"Unmapped transcript expression filter\n"
        f"{'=' * 50}\n"
        f"Input proteins:          {n_input:>10,}\n"
        f"Gene/intergenic (kept):  {n_gene_intergenic:>10,}\n"
        f"From unmapped tx:        {len(unmapped_proteins):>10,}\n"
        f"  Expressed (kept):      {n_kept_unmapped:>10,} "
        f"(TPM >= {args.min_tpm} in >= {args.min_samples} samples)\n"
        f"  Removed (no expr):     {n_removed:>10,}\n"
        f"{'─' * 50}\n"
        f"Output proteins:         {n_written:>10,}\n"
    )

    with open(args.stats, 'w') as sfh:
        sfh.write(stats)
    print(stats, file=sys.stderr)


if __name__ == '__main__':
    main()
