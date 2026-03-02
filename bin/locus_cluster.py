#!/usr/bin/env python3
"""
locus_cluster.py — Collapse transcripts that map to the same genomic locus.

Reads minimap2 PAF alignments, groups transcripts by overlapping genomic
coordinates (same chromosome, same strand), picks one representative per
locus (most aligned bases), and retains all unmapped transcripts.

Designed for cross-genus reference mapping where many transcripts won't
map, so unmapped transcripts are preserved by default.
"""

import argparse
import sys
from collections import defaultdict
from Bio import SeqIO


def parse_paf(paf_file, min_coverage, min_mapq):
    """Parse PAF file and return dict of query -> best alignment info.

    Filters on query coverage and mapping quality (mapq).  We do NOT
    filter on nucleotide identity (matches/block_len) because cross-genus
    divergence makes that metric unreliable — diverged orthologs can have
    <40% nt identity yet still map uniquely (mapq=60).

    PAF format (tab-separated):
        0  query name
        1  query length
        2  query start (0-based)
        3  query end
        4  strand (+/-)
        5  target name
        6  target length
        7  target start (0-based)
        8  target end
        9  number of matching bases
        10 alignment block length
        11 mapping quality (0-255)
    """
    alignments = {}

    with open(paf_file) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 12:
                continue

            query = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            target = fields[5]
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            block_len = int(fields[10])
            mapq = int(fields[11])

            # Coverage: fraction of query aligned
            coverage = (qend - qstart) / qlen if qlen > 0 else 0

            if coverage < min_coverage or mapq < min_mapq:
                continue

            # Keep best alignment per query (most matching bases)
            if query not in alignments or matches > alignments[query]['matches']:
                alignments[query] = {
                    'target': target,
                    'strand': strand,
                    'tstart': tstart,
                    'tend': tend,
                    'qlen': qlen,
                    'coverage': coverage,
                    'mapq': mapq,
                    'matches': matches,
                }

    return alignments


def merge_loci(alignments, max_gap):
    """Group overlapping/nearby alignments into genomic loci.

    Two alignments are in the same locus if they are on the same
    chromosome, same strand, and their genomic intervals overlap or
    are within max_gap bp of each other.
    """
    # Group by (target chromosome, strand)
    by_chrom_strand = defaultdict(list)
    for qname, aln in alignments.items():
        key = (aln['target'], aln['strand'])
        by_chrom_strand[key].append((aln['tstart'], aln['tend'], qname))

    loci = []
    locus_id = 0

    for (chrom, strand), intervals in sorted(by_chrom_strand.items()):
        intervals.sort()  # sort by start coordinate

        # Merge overlapping/nearby intervals
        merged_start = intervals[0][0]
        merged_end = intervals[0][1]
        current_members = [intervals[0][2]]

        for start, end, qname in intervals[1:]:
            if start <= merged_end + max_gap:
                # Overlaps or within gap tolerance — extend locus
                merged_end = max(merged_end, end)
                current_members.append(qname)
            else:
                # New locus
                locus_id += 1
                loci.append({
                    'locus_id': f'locus_{locus_id:06d}',
                    'chrom': chrom,
                    'strand': strand,
                    'start': merged_start,
                    'end': merged_end,
                    'members': current_members,
                })
                merged_start = start
                merged_end = end
                current_members = [qname]

        # Don't forget last locus in this chrom/strand
        locus_id += 1
        loci.append({
            'locus_id': f'locus_{locus_id:06d}',
            'chrom': chrom,
            'strand': strand,
            'start': merged_start,
            'end': merged_end,
            'members': current_members,
        })

    return loci


def pick_best_per_locus(loci, alignments):
    """Pick best transcript per locus (most aligned bases, then longest query).

    Returns:
        selected: dict of {transcript_id: locus_id} for representatives
        member_to_locus: dict of {transcript_id: locus_id} for all mapped
    """
    selected = {}
    member_to_locus = {}

    for locus in loci:
        members = locus['members']
        # Sort by: most matching bases, then longest query
        best = max(members, key=lambda q: (
            alignments[q]['matches'],
            alignments[q]['qlen'],
        ))
        selected[best] = locus['locus_id']
        for m in members:
            member_to_locus[m] = locus['locus_id']

    return selected, member_to_locus


def main():
    parser = argparse.ArgumentParser(
        description='Collapse transcripts by genomic locus (reference-guided)')
    parser.add_argument('--paf', required=True,
                        help='minimap2 PAF alignment file')
    parser.add_argument('--fasta', required=True,
                        help='Input representative FASTA')
    parser.add_argument('--out', required=True,
                        help='Output collapsed FASTA')
    parser.add_argument('--map', required=True,
                        help='Output locus map TSV')
    parser.add_argument('--stats', required=True,
                        help='Output stats file')
    parser.add_argument('--max-intron', type=int, default=200000,
                        help='Max gap between alignments for locus merging (bp) '
                             '[default: 200000]')
    parser.add_argument('--min-coverage', type=float, default=0.5,
                        help='Min query coverage to accept an alignment '
                             '[default: 0.5]')
    parser.add_argument('--min-mapq', type=int, default=5,
                        help='Min mapping quality (0-255) to accept an '
                             'alignment.  mapq=60 is unique; mapq>=5 '
                             'allows multi-mappers with a clear best '
                             '[default: 5]')
    args = parser.parse_args()

    # --- Parse PAF ---
    alignments = parse_paf(args.paf, args.min_coverage, args.min_mapq)
    print(f"Accepted alignments: {len(alignments)}", file=sys.stderr)

    # --- Read all transcript IDs ---
    all_ids = []
    seen = set()
    for rec in SeqIO.parse(args.fasta, 'fasta'):
        if rec.id not in seen:
            all_ids.append(rec.id)
            seen.add(rec.id)
    all_id_set = set(all_ids)
    mapped_ids = set(alignments.keys()) & all_id_set
    unmapped_ids = all_id_set - mapped_ids

    print(f"Total transcripts: {len(all_ids)}", file=sys.stderr)
    print(f"Mapped (accepted): {len(mapped_ids)}", file=sys.stderr)
    print(f"Unmapped/filtered: {len(unmapped_ids)}", file=sys.stderr)

    # --- Merge into loci ---
    # Only use alignments for transcripts in our FASTA
    filtered_alns = {k: v for k, v in alignments.items() if k in all_id_set}
    loci = merge_loci(filtered_alns, args.max_intron)

    # --- Pick best per locus ---
    selected, member_to_locus = pick_best_per_locus(loci, filtered_alns)

    # IDs to keep: best from each locus + all unmapped
    keep_ids = set(selected.keys()) | unmapped_ids

    print(f"Selected from loci: {len(selected)}", file=sys.stderr)
    print(f"Unmapped retained:  {len(unmapped_ids)}", file=sys.stderr)
    print(f"Total output:       {len(keep_ids)}", file=sys.stderr)

    # --- Write output FASTA (preserve input order) ---
    n_written = 0
    with open(args.out, 'w') as out_fh:
        for rec in SeqIO.parse(args.fasta, 'fasta'):
            if rec.id in keep_ids:
                SeqIO.write(rec, out_fh, 'fasta')
                n_written += 1

    # --- Write locus map ---
    with open(args.map, 'w') as map_fh:
        map_fh.write('transcript_id\tlocus_id\tstatus\n')
        for tid in all_ids:
            if tid in member_to_locus:
                status = 'selected' if tid in selected else 'collapsed'
                locus = member_to_locus[tid]
                map_fh.write(f'{tid}\t{locus}\t{status}\n')
            else:
                map_fh.write(f'{tid}\tunmapped\tretained\n')

    # --- Stats ---
    n_input = len(all_ids)
    n_mapped = len(mapped_ids)
    n_unmapped = len(unmapped_ids)
    n_loci = len(loci)
    n_multi = sum(1 for loc in loci if len(loc['members']) > 1)
    n_collapsed = n_input - n_written
    biggest_locus = max((len(loc['members']) for loc in loci), default=0)

    # Locus size distribution
    sizes = [len(loc['members']) for loc in loci]
    size_1 = sum(1 for s in sizes if s == 1)
    size_2_5 = sum(1 for s in sizes if 2 <= s <= 5)
    size_6_20 = sum(1 for s in sizes if 6 <= s <= 20)
    size_21plus = sum(1 for s in sizes if s > 20)

    stats = (
        f"Locus clustering (reference-guided)\n"
        f"{'=' * 42}\n"
        f"Input transcripts:     {n_input:>10,}\n"
        f"Mapped to genome:      {n_mapped:>10,} ({n_mapped/n_input*100:.1f}%)\n"
        f"Unmapped (retained):   {n_unmapped:>10,} ({n_unmapped/n_input*100:.1f}%)\n"
        f"{'─' * 42}\n"
        f"Genomic loci:          {n_loci:>10,}\n"
        f"  Singleton loci:      {size_1:>10,}\n"
        f"  2-5 transcripts:     {size_2_5:>10,}\n"
        f"  6-20 transcripts:    {size_6_20:>10,}\n"
        f"  >20 transcripts:     {size_21plus:>10,}\n"
        f"Largest locus:         {biggest_locus:>10,} transcripts\n"
        f"{'─' * 42}\n"
        f"Output transcripts:    {n_written:>10,}\n"
        f"Collapsed:             {n_collapsed:>10,} ({n_collapsed/n_input*100:.1f}%)\n"
    )

    with open(args.stats, 'w') as stats_fh:
        stats_fh.write(stats)

    print(stats, file=sys.stderr)


if __name__ == '__main__':
    main()
