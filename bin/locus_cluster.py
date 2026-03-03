#!/usr/bin/env python3
"""
locus_cluster.py — Collapse transcripts that map to the same genomic gene/locus.

Two modes:
  1. Gene-level (--gff provided): assign each mapped transcript to the
     reference gene it overlaps most, then pick one transcript per gene.
  2. Coordinate-overlap fallback (no --gff): merge overlapping alignment
     intervals into ad-hoc loci.

Unmapped transcripts are always retained.

Filters on mapping quality (mapq) and query coverage — NOT nucleotide
identity, which is unreliable for cross-genus mapping.
"""

import argparse
import bisect
import re
import sys
from collections import defaultdict
from Bio import SeqIO


# ── PAF parsing ──────────────────────────────────────────────────────

def parse_paf(paf_file, min_coverage, min_mapq):
    """Parse PAF, keep best hit per query.  Filter on coverage + mapq."""
    alignments = {}
    with open(paf_file) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 12:
                continue
            query  = f[0]
            qlen   = int(f[1])
            qstart = int(f[2])
            qend   = int(f[3])
            strand = f[4]
            target = f[5]
            tstart = int(f[7])
            tend   = int(f[8])
            matches = int(f[9])
            mapq   = int(f[11])

            coverage = (qend - qstart) / qlen if qlen > 0 else 0
            if coverage < min_coverage or mapq < min_mapq:
                continue

            if query not in alignments or matches > alignments[query]['matches']:
                alignments[query] = dict(
                    target=target, strand=strand, tstart=tstart, tend=tend,
                    qlen=qlen, coverage=coverage, mapq=mapq, matches=matches,
                )
    return alignments


# ── GFF gene parsing ────────────────────────────────────────────────

def parse_gff_genes(gff_file):
    """Parse gene records from GFF3, return sorted interval lists per
    (chrom, strand) keyed dict.

    Returns:
        genes_by_cs: dict of (chrom, strand) -> list of (start, end, gene_id)
                     sorted by start coordinate for bisect lookup.
    """
    genes_by_cs = defaultdict(list)
    id_re = re.compile(r'ID=([^;]+)')

    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'gene':
                continue
            chrom  = f[0]
            start  = int(f[3]) - 1   # GFF is 1-based → 0-based
            end    = int(f[4])        # GFF end is inclusive → exclusive
            strand = f[6]
            m = id_re.search(f[8])
            if not m:
                continue
            gene_id = m.group(1)
            # strip "gene-" prefix if present (EGAPx convention)
            if gene_id.startswith('gene-'):
                gene_id = gene_id[5:]
            genes_by_cs[(chrom, strand)].append((start, end, gene_id))

    # Sort by start for binary search
    for key in genes_by_cs:
        genes_by_cs[key].sort()

    n_genes = sum(len(v) for v in genes_by_cs.values())
    print(f"Parsed {n_genes} genes from GFF across "
          f"{len(genes_by_cs)} chrom/strand groups", file=sys.stderr)
    return genes_by_cs


def assign_to_genes(alignments, genes_by_cs):
    """Assign each alignment to the best-overlapping reference gene.

    For each mapped transcript, find all genes on the same chrom/strand
    that overlap the alignment interval, pick the gene with the largest
    overlap.  Transcripts that don't overlap any gene go to 'intergenic'.

    Returns:
        gene_members: dict of gene_id -> list of transcript IDs
        tx_to_gene:   dict of transcript_id -> gene_id
    """
    gene_members = defaultdict(list)
    tx_to_gene = {}
    n_intergenic = 0

    for qname, aln in alignments.items():
        key = (aln['target'], aln['strand'])
        genes = genes_by_cs.get(key, [])
        if not genes:
            # No genes on this chrom/strand — intergenic.
            # Each intergenic transcript gets its own group so it's retained.
            n_intergenic += 1
            igid = f'intergenic_{n_intergenic:06d}'
            gene_members[igid].append(qname)
            tx_to_gene[qname] = igid
            continue

        # Binary search: find genes whose start <= aln['tend']
        # and end >= aln['tstart'] (overlap condition)
        starts = [g[0] for g in genes]
        # Rightmost gene that could overlap: start < tend
        right_idx = bisect.bisect_left(starts, aln['tend'])

        best_gene = None
        best_overlap = 0

        # Scan backwards from right_idx
        for i in range(right_idx - 1, -1, -1):
            gstart, gend, gid = genes[i]
            if gend <= aln['tstart']:
                break  # no more overlap possible (sorted by start)
            # Compute overlap
            ovl_start = max(gstart, aln['tstart'])
            ovl_end = min(gend, aln['tend'])
            overlap = ovl_end - ovl_start
            if overlap > best_overlap:
                best_overlap = overlap
                best_gene = gid

        if best_gene:
            gene_members[best_gene].append(qname)
            tx_to_gene[qname] = best_gene
        else:
            # Mapped to a chrom with genes but didn't overlap any — intergenic
            n_intergenic += 1
            igid = f'intergenic_{n_intergenic:06d}'
            gene_members[igid].append(qname)
            tx_to_gene[qname] = igid

    return gene_members, tx_to_gene, n_intergenic


# ── Coordinate-overlap fallback ─────────────────────────────────────

def merge_loci(alignments, max_gap):
    """Group overlapping/nearby alignments into ad-hoc loci (no GFF)."""
    by_cs = defaultdict(list)
    for qname, aln in alignments.items():
        key = (aln['target'], aln['strand'])
        by_cs[key].append((aln['tstart'], aln['tend'], qname))

    gene_members = defaultdict(list)
    tx_to_gene = {}
    locus_id = 0

    for (chrom, strand), intervals in sorted(by_cs.items()):
        intervals.sort()
        merged_start = intervals[0][0]
        merged_end = intervals[0][1]
        current = [intervals[0][2]]

        for start, end, qname in intervals[1:]:
            if start <= merged_end + max_gap:
                merged_end = max(merged_end, end)
                current.append(qname)
            else:
                locus_id += 1
                lid = f'locus_{locus_id:06d}'
                gene_members[lid] = current
                for m in current:
                    tx_to_gene[m] = lid
                merged_start = start
                merged_end = end
                current = [qname]

        locus_id += 1
        lid = f'locus_{locus_id:06d}'
        gene_members[lid] = current
        for m in current:
            tx_to_gene[m] = lid

    return gene_members, tx_to_gene


# ── Best-per-group selection ────────────────────────────────────────

def pick_best_per_group(gene_members, alignments):
    """Pick best transcript per gene/locus: most aligned bases → longest.

    Returns:
        selected: dict of {transcript_id: group_id} for representatives
    """
    selected = {}
    for gid, members in gene_members.items():
        best = max(members, key=lambda q: (
            alignments[q]['matches'],
            alignments[q]['qlen'],
        ))
        selected[best] = gid
    return selected


# ── Main ────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Collapse transcripts by genomic gene/locus')
    parser.add_argument('--paf', required=True,
                        help='minimap2 PAF alignment file')
    parser.add_argument('--fasta', required=True,
                        help='Input representative FASTA')
    parser.add_argument('--out', required=True,
                        help='Output collapsed FASTA')
    parser.add_argument('--map', required=True,
                        help='Output map TSV')
    parser.add_argument('--stats', required=True,
                        help='Output stats file')
    parser.add_argument('--gff', default=None,
                        help='Reference GFF3 with gene annotations. '
                             'Enables gene-level collapse instead of '
                             'coordinate-overlap merging.')
    parser.add_argument('--max-intron', type=int, default=200000,
                        help='Max gap for coordinate-overlap merging (bp), '
                             'only used when --gff is not provided '
                             '[default: 200000]')
    parser.add_argument('--min-coverage', type=float, default=0.5,
                        help='Min query coverage [default: 0.5]')
    parser.add_argument('--min-mapq', type=int, default=5,
                        help='Min mapping quality [default: 5]')
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
    filtered_alns = {k: v for k, v in alignments.items() if k in all_id_set}

    print(f"Total transcripts: {len(all_ids)}", file=sys.stderr)
    print(f"Mapped (accepted): {len(mapped_ids)}", file=sys.stderr)
    print(f"Unmapped/filtered: {len(unmapped_ids)}", file=sys.stderr)

    # --- Group transcripts by gene/locus ---
    if args.gff:
        genes_by_cs = parse_gff_genes(args.gff)
        gene_members, tx_to_gene, n_intergenic = assign_to_genes(filtered_alns, genes_by_cs)
        mode = 'gene'
    else:
        gene_members, tx_to_gene = merge_loci(filtered_alns, args.max_intron)
        mode = 'locus'
        n_intergenic = 0

    # --- Pick best per group ---
    selected = pick_best_per_group(gene_members, filtered_alns)

    # IDs to keep: best from each gene/locus + all unmapped
    keep_ids = set(selected.keys()) | unmapped_ids

    # Gene groups = total groups minus intergenic singletons
    n_gene_groups = len(gene_members) - n_intergenic

    print(f"Groups ({mode}):       {n_gene_groups}", file=sys.stderr)
    if mode == 'gene':
        print(f"Intergenic:          {n_intergenic}", file=sys.stderr)
    print(f"Selected:            {len(selected)}", file=sys.stderr)
    print(f"Unmapped retained:   {len(unmapped_ids)}", file=sys.stderr)
    print(f"Total output:        {len(keep_ids)}", file=sys.stderr)

    # --- Write output FASTA ---
    n_written = 0
    with open(args.out, 'w') as out_fh:
        for rec in SeqIO.parse(args.fasta, 'fasta'):
            if rec.id in keep_ids:
                SeqIO.write(rec, out_fh, 'fasta')
                n_written += 1

    # --- Write map ---
    with open(args.map, 'w') as map_fh:
        map_fh.write(f'transcript_id\t{mode}_id\tstatus\n')
        for tid in all_ids:
            if tid in tx_to_gene:
                status = 'selected' if tid in selected else 'collapsed'
                map_fh.write(f'{tid}\t{tx_to_gene[tid]}\t{status}\n')
            else:
                map_fh.write(f'{tid}\tunmapped\tretained\n')

    # --- Stats ---
    n_input = len(all_ids)
    n_mapped = len(mapped_ids)
    n_unmapped = len(unmapped_ids)
    n_collapsed = n_input - n_written

    sizes = [len(v) for k, v in gene_members.items()
             if not k.startswith('intergenic_')]
    size_1 = sum(1 for s in sizes if s == 1)
    size_2_5 = sum(1 for s in sizes if 2 <= s <= 5)
    size_6_20 = sum(1 for s in sizes if 6 <= s <= 20)
    size_21plus = sum(1 for s in sizes if s > 20)
    biggest = max(sizes, default=0)

    header = f"Gene-level collapse (GFF)" if mode == 'gene' else \
             f"Locus clustering (coordinate overlap)"

    stats = (
        f"{header}\n"
        f"{'=' * 50}\n"
        f"Input transcripts:       {n_input:>10,}\n"
        f"Mapped to genome:        {n_mapped:>10,} "
        f"({n_mapped/n_input*100:.1f}%)\n"
        f"Unmapped (retained):     {n_unmapped:>10,} "
        f"({n_unmapped/n_input*100:.1f}%)\n"
    )
    if mode == 'gene':
        stats += (
            f"Assigned to genes:       {n_mapped - n_intergenic:>10,}\n"
            f"Intergenic (no gene):    {n_intergenic:>10,}\n"
        )
    stats += (
        f"{'─' * 50}\n"
        f"Reference genes hit:     {n_gene_groups:>10,}\n"
        f"  1 transcript:          {size_1:>10,}\n"
        f"  2-5 transcripts:       {size_2_5:>10,}\n"
        f"  6-20 transcripts:      {size_6_20:>10,}\n"
        f"  >20 transcripts:       {size_21plus:>10,}\n"
        f"Largest gene group:      {biggest:>10,} transcripts\n"
        f"{'─' * 50}\n"
        f"Output transcripts:      {n_written:>10,}\n"
        f"Collapsed:               {n_collapsed:>10,} "
        f"({n_collapsed/n_input*100:.1f}%)\n"
    )

    with open(args.stats, 'w') as stats_fh:
        stats_fh.write(stats)

    print(stats, file=sys.stderr)


if __name__ == '__main__':
    main()
