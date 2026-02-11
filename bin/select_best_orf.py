#!/usr/bin/env python3
"""Select the single best ORF per SuperTranscript gene.

Usage:
    select_best_orf.py <psauron_score.csv> <TD2.pep> <TD2.gff3> <species_label>

Outputs:
    - {species_label}.faa    : proteins with gene-level IDs
    - best_orfs.gff3         : filtered GFF3 for retained ORFs only
    - orf_to_gene_map.tsv    : gene_id -> orf_id -> psauron_score -> orf_length
"""

import csv
import sys
from collections import defaultdict
from Bio import SeqIO

psauron_csv = sys.argv[1]
td2_pep     = sys.argv[2]
td2_gff3    = sys.argv[3]
species     = sys.argv[4]

# --- Parse PSAURON scores ---
scores = {}
with open(psauron_csv) as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    for row in reader:
        scores[row[0]] = float(row[1])

# --- Group ORFs by parent SuperTranscript ---
orfs_by_gene = defaultdict(list)
for rec in SeqIO.parse(td2_pep, "fasta"):
    orf_id = rec.id
    # Parent gene = everything before the last .pN
    gene_id = orf_id.rsplit(".", 1)[0]
    psauron = scores.get(orf_id, 0.0)
    orfs_by_gene[gene_id].append((psauron, len(rec.seq), orf_id, rec))

# --- Select best: highest PSAURON, then longest ---
kept_orfs = set()
with open("orf_to_gene_map.tsv", "w") as mapf:
    mapf.write("gene_id\torf_id\tpsauron_score\torf_length\n")
    with open(f"{species}.faa", "w") as out:
        for gene_id in sorted(orfs_by_gene):
            orfs = orfs_by_gene[gene_id]
            orfs.sort(key=lambda x: (x[0], x[1]), reverse=True)
            best_psauron, best_len, best_orf_id, best_rec = orfs[0]
            out.write(f">{gene_id}\n{str(best_rec.seq)}\n")
            mapf.write(f"{gene_id}\t{best_orf_id}\t{best_psauron}\t{best_len}\n")
            kept_orfs.add(best_orf_id)

# --- Filter GFF3 ---
with open(td2_gff3) as gff_in, open("best_orfs.gff3", "w") as gff_out:
    for line in gff_in:
        if line.startswith("#"):
            gff_out.write(line)
            continue
        if any(orf_id in line for orf_id in kept_orfs):
            gff_out.write(line)

n_genes = len(orfs_by_gene)
n_total = sum(len(v) for v in orfs_by_gene.values())
print(f"Selected {n_genes} best proteins from {n_total} total ORFs "
      f"({n_total - n_genes} secondary ORFs removed)")
