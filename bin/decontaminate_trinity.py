#!/usr/bin/env python3
"""Extract original Trinity transcripts belonging to taxonomy-clean clusters.

Usage:
    decontaminate_trinity.py <trinity.fasta> <corset_clusters.txt> \
        <filtered_reps.fasta> <output.fasta> <output_stats.txt>

Inputs:
    trinity.fasta          — original Trinity assembly (all transcripts)
    corset_clusters.txt    — Corset cluster assignments (transcript → cluster)
    filtered_reps.fasta    — taxonomy-filtered representative FASTA from
                             MMSEQS2_TAXONOMY (headers = cluster IDs)

Logic:
    1. Read cluster IDs from filtered_reps.fasta headers (= clean set)
    2. Read Corset mapping to find all Trinity transcripts in clean clusters
    3. Write those transcripts from the original Trinity FASTA

This produces a decontaminated Trinity assembly: all original isoforms
from clusters classified as Viridiplantae or no-hit, preserving the
full isoform diversity of the original assembly minus contaminants.
"""

import sys
from Bio import SeqIO

if len(sys.argv) != 6:
    print(f"Usage: {sys.argv[0]} <trinity.fasta> <corset_clusters.txt> "
          f"<filtered_reps.fasta> <output.fasta> <output_stats.txt>",
          file=sys.stderr)
    sys.exit(1)

trinity_fasta   = sys.argv[1]
clusters_file   = sys.argv[2]
filtered_fasta  = sys.argv[3]
out_fasta       = sys.argv[4]
out_stats       = sys.argv[5]

# 1. Collect clean cluster IDs from taxonomy-filtered representatives
clean_clusters = set()
for rec in SeqIO.parse(filtered_fasta, "fasta"):
    clean_clusters.add(rec.id)

print(f"Clean clusters (passed taxonomy filter): {len(clean_clusters):,}")

# 2. Read Corset cluster assignments and identify clean transcripts
trans_to_cluster = {}
with open(clusters_file) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            trans_to_cluster[parts[0]] = parts[1]

clean_transcripts = set()
contam_transcripts = set()
for tid, cid in trans_to_cluster.items():
    if cid in clean_clusters:
        clean_transcripts.add(tid)
    else:
        contam_transcripts.add(tid)

print(f"Transcripts in clean clusters: {len(clean_transcripts):,}")
print(f"Transcripts in contaminant clusters: {len(contam_transcripts):,}")

# 3. Extract clean transcripts from original Trinity FASTA
n_written = 0
n_skipped = 0
n_not_clustered = 0

with open(out_fasta, "w") as fout:
    for rec in SeqIO.parse(trinity_fasta, "fasta"):
        if rec.id in clean_transcripts:
            fout.write(f">{rec.id} {rec.description.split(' ', 1)[1] if ' ' in rec.description else ''}\n")
            fout.write(f"{str(rec.seq)}\n")
            n_written += 1
        elif rec.id in contam_transcripts:
            n_skipped += 1
        else:
            # Transcript not in any Corset cluster (should be rare)
            n_not_clustered += 1

total = n_written + n_skipped + n_not_clustered
pct_kept = n_written / total * 100 if total else 0
pct_removed = n_skipped / total * 100 if total else 0

print(f"\nDecontaminated Trinity assembly:")
print(f"  Kept:           {n_written:,} ({pct_kept:.1f}%)")
print(f"  Removed:        {n_skipped:,} ({pct_removed:.1f}%)")
if n_not_clustered:
    print(f"  Not clustered:  {n_not_clustered:,} (excluded)")
print(f"  Total input:    {total:,}")

# 4. Write stats file
with open(out_stats, "w") as sf:
    sf.write("Decontaminated Trinity Assembly Statistics\n")
    sf.write("=" * 45 + "\n\n")
    sf.write(f"Clean clusters (Viridiplantae + no-hit): {len(clean_clusters):,}\n")
    sf.write(f"Contaminant clusters (removed):          "
             f"{len(trans_to_cluster) - len(clean_transcripts) - len(contam_transcripts) + len(contam_transcripts):,}\n\n")
    n_contam_clusters = len(set(trans_to_cluster.values()) - clean_clusters)
    sf.write(f"Contaminant clusters:  {n_contam_clusters:,}\n\n")
    sf.write(f"Trinity transcripts kept:     {n_written:>10,} ({pct_kept:.1f}%)\n")
    sf.write(f"Trinity transcripts removed:  {n_skipped:>10,} ({pct_removed:.1f}%)\n")
    if n_not_clustered:
        sf.write(f"Trinity transcripts unclustered: {n_not_clustered:>7,}\n")
    sf.write(f"Total input transcripts:      {total:>10,}\n")
