#!/usr/bin/env python3
"""Select best MetaEuk protein per gene from easy-predict output.

Usage:
    metaeuk_select_best.py <metaeuk.fas> <output.faa> <output_map.tsv>

MetaEuk easy-predict .fas header format:
    >swissprot_id|gene_id|strand|score|evalue|num_exons

Selection: highest score -> lowest evalue -> longest protein.
"""

import sys
from collections import defaultdict

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <metaeuk.fas> <output.faa> <output_map.tsv>",
          file=sys.stderr)
    sys.exit(1)

infile, out_faa, out_map = sys.argv[1], sys.argv[2], sys.argv[3]

genes = defaultdict(list)
with open(infile) as f:
    header = None
    seq_lines = []
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                seq = "".join(seq_lines)
                genes[gene_id].append((score, evalue, len(seq), prot_id, seq))
            header = line[1:]
            seq_lines = []
            parts = header.split("|")
            swissprot_id = parts[0]
            gene_id = parts[1]
            score = float(parts[3])
            evalue = float(parts[4])
            prot_id = f"{gene_id}_{swissprot_id}"
        else:
            seq_lines.append(line)
    if header:
        seq = "".join(seq_lines)
        genes[gene_id].append((score, evalue, len(seq), prot_id, seq))

n_total = sum(len(v) for v in genes.values())

with open(out_faa, "w") as faa, open(out_map, "w") as mapf:
    mapf.write("gene_id\tprot_id\tscore\tevalue\tprot_length\n")
    for gene_id in sorted(genes):
        candidates = genes[gene_id]
        candidates.sort(key=lambda x: (x[0], -x[1], x[2]), reverse=True)
        best_score, best_eval, best_len, best_prot_id, best_seq = candidates[0]
        faa.write(f">{gene_id}\n{best_seq}\n")
        mapf.write(f"{gene_id}\t{best_prot_id}\t{best_score}\t{best_eval}\t{best_len}\n")

print(f"MetaEuk: {len(genes)} best proteins from {n_total} total predictions "
      f"({n_total - len(genes)} secondary removed)")
