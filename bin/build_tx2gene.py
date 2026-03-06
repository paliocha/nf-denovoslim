#!/usr/bin/env python3
"""Build a tx2gene mapping for tximport gene-level count aggregation.

Chains the pipeline's collapsing steps to map every Trinity transcript
in the decontaminated assembly to its final gene ID:

    Trinity transcript  →  Corset cluster  →  nt dedup representative
                        →  locus representative  →  protein dedup representative

Usage:
    build_tx2gene.py \\
        --corset-clusters clusters.txt \\
        --decontam-fasta decontaminated.fasta \\
        --nt-dedup-tsv clust_cluster.tsv \\
        [--locus-map locus_map.tsv] \\
        [--protein-dedup-map protein_cluster_map.tsv] \\
        --output tx2gene.tsv \\
        --stats tx2gene_stats.txt

Outputs:
    tx2gene.tsv  — two-column (TXNAME, GENEID) for tximport::tximport()
    tx2gene_stats.txt — summary of mapping chain resolution
"""

import argparse
import os
import sys
from collections import defaultdict


def count_fasta_ids(path):
    """Return set of sequence IDs from a FASTA file."""
    ids = set()
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                ids.add(line[1:].split()[0])
    return ids


def parse_corset_clusters(path):
    """Parse Corset cluster file: transcript_id → cluster_id."""
    tx_to_cluster = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                tx_to_cluster[parts[0]] = parts[1]
    return tx_to_cluster


def parse_nt_dedup_tsv(path):
    """Parse MMseqs2 easy-cluster TSV: build member → representative map.

    MMseqs2 clust_cluster.tsv format: representative_id<TAB>member_id
    One row per member (including self-hits for representatives).
    """
    member_to_rep = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                rep_id = parts[0].split()[0]  # strip trailing description
                mem_id = parts[1].split()[0]
                member_to_rep[mem_id] = rep_id
    return member_to_rep


def parse_locus_map(path):
    """Parse locus_cluster.py map: build collapsed → selected representative.

    Format: transcript_id<TAB>gene_or_locus_id<TAB>status
    Status values: 'selected', 'collapsed', 'retained'

    For each group, find the 'selected' member; map all 'collapsed' members
    to that selected ID.
    """
    # First pass: collect groups and find selected member per group
    group_members = defaultdict(list)  # group → [(tid, status)]
    with open(path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                tid, group_id, status = parts[0], parts[1], parts[2]
                group_members[group_id].append((tid, status))

    # Build mapping: collapsed → selected
    collapsed_to_selected = {}
    for group_id, members in group_members.items():
        selected = [tid for tid, status in members if status == 'selected']
        if not selected:
            continue  # unmapped / retained — keep own ID
        rep = selected[0]
        for tid, status in members:
            if status == 'collapsed':
                collapsed_to_selected[tid] = rep

    return collapsed_to_selected


def parse_protein_dedup_map(path):
    """Parse protein dedup cluster map: gene_id → protein_group_id.

    Format (with header): gene_id<TAB>protein_group_id
    """
    gene_to_group = {}
    with open(path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gene_to_group[parts[0]] = parts[1]
    return gene_to_group


def main():
    parser = argparse.ArgumentParser(
        description='Build tx2gene mapping for tximport')
    parser.add_argument('--corset-clusters', required=True,
                        help='Corset cluster file (transcript → cluster)')
    parser.add_argument('--decontam-fasta', required=True,
                        help='Decontaminated Trinity FASTA')
    parser.add_argument('--nt-dedup-tsv', required=True,
                        help='MMseqs2 nt dedup cluster TSV (rep → member)')
    parser.add_argument('--locus-map', default='',
                        help='Locus cluster map TSV (optional)')
    parser.add_argument('--protein-dedup-map', default='',
                        help='Protein dedup cluster map TSV (optional)')
    parser.add_argument('--output', required=True,
                        help='Output tx2gene.tsv')
    parser.add_argument('--stats', required=True,
                        help='Output stats file')
    args = parser.parse_args()

    # --- Load mappings ---
    print("Loading Corset clusters...")
    tx_to_cluster = parse_corset_clusters(args.corset_clusters)
    print(f"  {len(tx_to_cluster):,} transcripts in {len(set(tx_to_cluster.values())):,} clusters")

    print("Loading nt dedup cluster map...")
    nt_member_to_rep = parse_nt_dedup_tsv(args.nt_dedup_tsv)
    n_nt_reps = len(set(nt_member_to_rep.values()))
    print(f"  {len(nt_member_to_rep):,} members → {n_nt_reps:,} representatives")

    locus_collapsed_to_selected = {}
    if args.locus_map and os.path.isfile(args.locus_map):
        print("Loading locus cluster map...")
        locus_collapsed_to_selected = parse_locus_map(args.locus_map)
        print(f"  {len(locus_collapsed_to_selected):,} collapsed → selected mappings")

    protein_gene_to_group = {}
    if args.protein_dedup_map and os.path.isfile(args.protein_dedup_map):
        print("Loading protein dedup map...")
        protein_gene_to_group = parse_protein_dedup_map(args.protein_dedup_map)
        n_prot_groups = len(set(protein_gene_to_group.values()))
        print(f"  {len(protein_gene_to_group):,} genes → {n_prot_groups:,} protein groups")

    # --- Get decontaminated transcript IDs ---
    print("Reading decontaminated Trinity FASTA IDs...")
    decontam_ids = count_fasta_ids(args.decontam_fasta)
    print(f"  {len(decontam_ids):,} transcripts")

    # --- Resolve chain for each transcript ---
    print("Resolving tx2gene chain...")

    n_no_cluster = 0
    n_total = 0
    gene_ids = set()

    with open(args.output, 'w') as out:
        out.write("TXNAME\tGENEID\n")

        for tx_id in sorted(decontam_ids):
            # Step 1: transcript → Corset cluster
            cluster = tx_to_cluster.get(tx_id)
            if cluster is None:
                n_no_cluster += 1
                continue

            gene = cluster

            # Step 2: Corset cluster → nt dedup representative
            gene = nt_member_to_rep.get(gene, gene)

            # Step 3: nt dedup rep → locus cluster representative
            if locus_collapsed_to_selected:
                gene = locus_collapsed_to_selected.get(gene, gene)

            # Step 4: locus rep → protein dedup representative
            # Note: only genes WITH proteins appear in protein_dedup_map.
            # Genes without proteins keep their locus/nt-dedup gene ID.
            if protein_gene_to_group:
                gene = protein_gene_to_group.get(gene, gene)

            out.write(f"{tx_id}\t{gene}\n")
            gene_ids.add(gene)
            n_total += 1

    # --- Stats ---
    with open(args.stats, 'w') as stats:
        stats.write("tx2gene mapping statistics\n")
        stats.write("=" * 50 + "\n")
        stats.write(f"Decontaminated transcripts:  {len(decontam_ids):,}\n")
        stats.write(f"Mapped to genes:             {n_total:,}\n")
        stats.write(f"No Corset cluster:           {n_no_cluster:,}\n")
        stats.write(f"Unique gene IDs:             {len(gene_ids):,}\n")
        stats.write(f"\nChain resolution steps:\n")
        stats.write(f"  Corset clusters:           {len(set(tx_to_cluster.values())):,}\n")
        stats.write(f"  After nt dedup:            {n_nt_reps:,}\n")
        if locus_collapsed_to_selected:
            n_after_locus = n_nt_reps - len(locus_collapsed_to_selected)
            stats.write(f"  After locus clustering:    ~{n_after_locus:,}\n")
        if protein_gene_to_group:
            stats.write(f"  After protein dedup:       {len(set(protein_gene_to_group.values())):,}\n")
        stats.write(f"\nMean transcripts per gene:   {n_total / len(gene_ids):.1f}\n"
                     if gene_ids else "")

    print(f"tx2gene written: {n_total:,} transcripts → {len(gene_ids):,} genes")
    if n_no_cluster:
        print(f"WARNING: {n_no_cluster:,} transcripts had no Corset cluster", file=sys.stderr)


if __name__ == "__main__":
    main()
