#!/usr/bin/env python3
"""
miniprot_filter.py — Filter proteins based on miniprot genome validation.

Reads a miniprot GFF and the input protein FASTA.  Keeps proteins that:
  1. Map to the genome with Identity >= --min-identity, OR
  2. Are >= --min-length amino acids (trusted by size alone)

Short proteins that lack genomic support are removed as potential
prediction artifacts.
"""

import argparse
import sys
from Bio import SeqIO


def parse_miniprot_gff(gff_path):
    """Extract best alignment identity per protein from miniprot GFF.

    Parses mRNA lines for Target= (protein ID) and Identity= (fraction).
    Returns dict mapping protein_id -> best_identity.
    """
    best = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'mRNA':
                continue
            target = None
            identity = 0.0
            for attr in fields[8].split(';'):
                if attr.startswith('Target='):
                    # Target=PROT_ID start end [+-]
                    target = attr.split('=', 1)[1].split()[0]
                elif attr.startswith('Identity='):
                    identity = float(attr.split('=', 1)[1])
            if target is not None:
                if target not in best or identity > best[target]:
                    best[target] = identity
    return best


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--proteins', required=True,
                        help='Input protein FASTA')
    parser.add_argument('--gff', required=True,
                        help='Miniprot GFF output')
    parser.add_argument('--min-identity', type=float, default=0.5,
                        help='Min alignment identity to count as '
                             'genome-validated (default: 0.5)')
    parser.add_argument('--min-length', type=int, default=75,
                        help='Min protein length (aa) to keep without '
                             'genome validation (default: 75)')
    parser.add_argument('--out', required=True,
                        help='Output filtered FASTA')
    parser.add_argument('--stats', required=True,
                        help='Output stats file')
    args = parser.parse_args()

    # Parse GFF for best identity per protein
    best_identity = parse_miniprot_gff(args.gff)

    # Filter proteins
    n_total = 0
    n_kept = 0
    n_kept_by_both = 0       # long AND mapped
    n_kept_by_length = 0     # long but unmapped / low-identity
    n_kept_by_mapping = 0    # short but genome-validated
    n_dropped = 0

    kept_records = []
    for record in SeqIO.parse(args.proteins, 'fasta'):
        n_total += 1
        prot_len = len(record.seq)
        is_long = prot_len >= args.min_length
        is_validated = (record.id in best_identity
                        and best_identity[record.id] >= args.min_identity)

        if is_long and is_validated:
            n_kept_by_both += 1
            kept_records.append(record)
            n_kept += 1
        elif is_long:
            n_kept_by_length += 1
            kept_records.append(record)
            n_kept += 1
        elif is_validated:
            n_kept_by_mapping += 1
            kept_records.append(record)
            n_kept += 1
        else:
            n_dropped += 1

    # Write filtered FASTA
    SeqIO.write(kept_records, args.out, 'fasta')

    # Write stats
    n_mapped = len(best_identity)
    n_mapped_pass = sum(1 for v in best_identity.values()
                        if v >= args.min_identity)
    pct_dropped = n_dropped / n_total * 100 if n_total > 0 else 0

    with open(args.stats, 'w') as fh:
        fh.write(f'Miniprot genome validation '
                 f'(identity >= {args.min_identity}, '
                 f'length >= {args.min_length} aa)\n')
        fh.write(f'Input proteins:             {n_total}\n')
        fh.write(f'Mapped to genome:           {n_mapped}\n')
        fh.write(f'Mapped >= identity cutoff:  {n_mapped_pass}\n')
        fh.write(f'Kept (total):               {n_kept}\n')
        fh.write(f'  Kept (long + mapped):     {n_kept_by_both}\n')
        fh.write(f'  Kept (long only):         {n_kept_by_length}\n')
        fh.write(f'  Kept (mapped only):       {n_kept_by_mapping}\n')
        fh.write(f'Dropped:                    {n_dropped} ({pct_dropped:.1f}%)\n')

    print(f'Miniprot filter: {n_total} -> {n_kept} '
          f'({n_dropped} dropped, {pct_dropped:.1f}%)')


if __name__ == '__main__':
    main()
