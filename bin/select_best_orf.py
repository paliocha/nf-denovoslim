#!/usr/bin/env python3
"""Select the single best ORF per representative gene.

Usage:
    select_best_orf.py <psauron_score.csv> <TD2.pep> <TD2.gff3> <species_label>

Outputs:
    - {species_label}.faa    : proteins with gene-level IDs
    - best_orfs.gff3         : filtered GFF3 for retained ORFs only
    - orf_to_gene_map.tsv    : gene_id -> orf_id -> psauron_score -> orf_length -> completeness
"""

import csv
import re
import sys
from collections import defaultdict
from Bio import SeqIO

if len(sys.argv) != 5:
    print(f"Usage: {sys.argv[0]} <psauron_score.csv> <TD2.pep> <TD2.gff3> <species_label>",
          file=sys.stderr)
    sys.exit(1)

psauron_csv = sys.argv[1]
td2_pep     = sys.argv[2]
td2_gff3    = sys.argv[3]
species     = sys.argv[4]

# --- Parse PSAURON scores ---
# PSAURON writes preamble lines (command echo, summary stats) before the CSV.
# Scan for the header line, then use DictReader for column-name stability.
scores = {}
with open(psauron_csv) as f:
    # Skip preamble until we find the CSV header
    for line in f:
        if line.startswith("description,"):
            break
    else:
        print("ERROR: could not find 'description,' header in PSAURON CSV", file=sys.stderr)
        sys.exit(1)
    reader = csv.DictReader(f, fieldnames=line.strip().split(","))
    for row in reader:
        orf_id = row["description"]
        try:
            scores[orf_id] = float(row["in-frame_score"])
        except (ValueError, KeyError):
            # Skip malformed rows (e.g. trailing blank lines)
            continue

# --- Parse completeness from TD2 .pep headers ---
# TD2 .pep description contains "type:complete", "type:5prime_partial",
# "type:3prime_partial", or "type:internal".
TYPE_RE = re.compile(r'type:(\S+)')

def parse_completeness(description):
    """Extract ORF type from TD2 .pep description."""
    m = TYPE_RE.search(description)
    return m.group(1) if m else "unknown"

# --- Group ORFs by parent gene ---
orfs_by_gene = defaultdict(list)
for rec in SeqIO.parse(td2_pep, "fasta"):
    orf_id = rec.id
    # Parent gene = everything before the last .pN
    gene_id = orf_id.rsplit(".", 1)[0]
    psauron = scores.get(orf_id, 0.0)
    completeness = parse_completeness(rec.description)
    orfs_by_gene[gene_id].append((psauron, len(rec.seq), orf_id, rec, completeness))

# --- Select best: highest PSAURON, then longest ---
kept_orfs = set()
with open("orf_to_gene_map.tsv", "w") as mapf:
    mapf.write("gene_id\torf_id\tpsauron_score\torf_length\tcompleteness\n")
    with open(f"{species}.faa", "w") as out:
        for gene_id in sorted(orfs_by_gene):
            orfs = orfs_by_gene[gene_id]
            orfs.sort(key=lambda x: (x[0], x[1]), reverse=True)
            best_psauron, best_len, best_orf_id, best_rec, best_compl = orfs[0]
            out.write(f">{gene_id}\n{str(best_rec.seq)}\n")
            mapf.write(f"{gene_id}\t{best_orf_id}\t{best_psauron}\t{best_len}\t{best_compl}\n")
            kept_orfs.add(best_orf_id)

# --- Filter GFF3 ---
# Parse the ID= attribute to avoid substring false positives
# (e.g. gene1.p1 matching gene1.p10)
id_pattern = re.compile(r'ID=([^;]+)')
parent_pattern = re.compile(r'Parent=([^;]+)')
with open(td2_gff3) as gff_in, open("best_orfs.gff3", "w") as gff_out:
    for line in gff_in:
        if line.startswith("#"):
            gff_out.write(line)
            continue
        id_match = id_pattern.search(line)
        parent_match = parent_pattern.search(line)
        gff_id = id_match.group(1) if id_match else None
        gff_parent = parent_match.group(1) if parent_match else None
        if (gff_id and gff_id in kept_orfs) or (gff_parent and gff_parent in kept_orfs):
            gff_out.write(line)

n_genes = len(orfs_by_gene)
n_total = sum(len(v) for v in orfs_by_gene.values())
print(f"Selected {n_genes} best proteins from {n_total} total ORFs "
      f"({n_total - n_genes} secondary ORFs removed)")
