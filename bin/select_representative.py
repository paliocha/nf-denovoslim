#!/usr/bin/env python3
"""Select the longest Trinity transcript per Corset cluster.

Usage:
    select_representative.py <trinity.fasta> <corset_clusters.txt> <output.fasta> <output_map.tsv>

Corset clusters format (tab-separated):
    transcript_id    cluster_id

Outputs:
    - output.fasta   : one representative sequence per cluster (longest wins),
                        header is the cluster_id
    - output_map.tsv : cluster_id -> representative_transcript_id -> length
"""

import sys
from collections import defaultdict
from Bio import SeqIO

if len(sys.argv) != 5:
    print(f"Usage: {sys.argv[0]} <trinity.fasta> <corset_clusters.txt> "
          f"<output.fasta> <output_map.tsv>", file=sys.stderr)
    sys.exit(1)

trinity_fasta = sys.argv[1]
clusters_file = sys.argv[2]
out_fasta     = sys.argv[3]
out_map       = sys.argv[4]

# --- Parse Corset cluster assignments ---
# transcript_id → cluster_id
trans_to_cluster = {}
with open(clusters_file) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            trans_to_cluster[parts[0]] = parts[1]

clusters = defaultdict(list)
for tid, cid in trans_to_cluster.items():
    clusters[cid].append(tid)

print(f"Parsed {len(trans_to_cluster):,} transcripts in {len(clusters):,} clusters")

# --- Index Trinity FASTA ---
# Build a dict of transcript_id → (length, record) for lookup
seq_index = {}
for rec in SeqIO.parse(trinity_fasta, "fasta"):
    seq_index[rec.id] = rec

print(f"Indexed {len(seq_index):,} sequences from Trinity FASTA")

# --- Select longest per cluster ---
n_singletons = 0
n_multi = 0
n_missing = 0

with open(out_fasta, "w") as fasta_out, open(out_map, "w") as map_out:
    map_out.write("cluster_id\trepresentative_transcript\tlength\tn_transcripts\n")

    for cid in sorted(clusters):
        transcripts = clusters[cid]
        best_rec = None
        best_len = -1

        for tid in transcripts:
            if tid in seq_index:
                rec = seq_index[tid]
                if len(rec.seq) > best_len:
                    best_len = len(rec.seq)
                    best_rec = rec

        if best_rec is None:
            n_missing += 1
            continue

        # Write with cluster ID as header
        fasta_out.write(f">{cid}\n{str(best_rec.seq)}\n")
        map_out.write(f"{cid}\t{best_rec.id}\t{best_len}\t{len(transcripts)}\n")

        if len(transcripts) == 1:
            n_singletons += 1
        else:
            n_multi += 1

total = n_singletons + n_multi
print(f"Selected {total:,} representatives "
      f"({n_singletons:,} singletons, {n_multi:,} multi-transcript clusters)")
if n_missing > 0:
    print(f"WARNING: {n_missing:,} clusters had no matching Trinity sequences",
          file=sys.stderr)
