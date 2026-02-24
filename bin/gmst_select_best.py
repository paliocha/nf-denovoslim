#!/usr/bin/env python3
"""Select best GeneMarkS-T ORF per gene.

Usage:
    gmst_select_best.py <input.fasta> <gmst.faa> <gmst.fnn> \
        <output.faa> <output.fnn> <output_map.tsv>

Reads the set of known gene IDs from the input FASTA, then maps
GeneMarkS-T ORF headers back to source genes.  Selects one best ORF
per gene: complete > partial, then longest protein.
"""

import sys
from collections import defaultdict

if len(sys.argv) != 7:
    print(f"Usage: {sys.argv[0]} <input.fasta> <gmst.faa> <gmst.fnn> "
          f"<output.faa> <output.fnn> <output_map.tsv>", file=sys.stderr)
    sys.exit(1)

input_fasta, gmst_faa, gmst_fnn = sys.argv[1], sys.argv[2], sys.argv[3]
out_faa, out_fnn, out_map = sys.argv[4], sys.argv[5], sys.argv[6]


# Build set of known gene IDs from input FASTA
known_genes = set()
with open(input_fasta) as f:
    for line in f:
        if line.startswith(">"):
            known_genes.add(line[1:].strip().split()[0])


def read_fasta(path):
    """Read FASTA file, return list of (header, sequence) tuples."""
    seqs = []
    hdr = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr is not None:
                    seqs.append((hdr, "".join(parts)))
                hdr = line[1:]
                parts = []
            else:
                parts.append(line)
        if hdr is not None:
            seqs.append((hdr, "".join(parts)))
    return seqs


def extract_gene_id(header, known):
    """Extract the source gene ID from a GeneMarkS-T ORF header.

    Tries multiple strategies:
      1. First token with trailing _N stripped (most gmst.pl versions)
      2. '|'-delimited: field containing a known gene ID
      3. Longest known ID that is a substring of the token (fallback)
    """
    token = header.split()[0]

    # Strategy 1: strip trailing _digits
    parts = token.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit() and parts[0] in known:
        return parts[0]

    # Strategy 2: pipe-delimited fields
    for field in token.split("|"):
        if field in known:
            return field

    # Strategy 3: longest known ID that is a substring of the token
    best = None
    for gid in known:
        if gid in token:
            if best is None or len(gid) > len(best):
                best = gid
    if best:
        return best

    # Last resort: use full first token (may fail downstream)
    return token


def infer_completeness(header, seq):
    """Infer ORF completeness from header annotations or sequence."""
    h_lower = header.lower()
    if "complete" in h_lower and "partial" not in h_lower:
        return "complete"
    if "partial" in h_lower:
        for tag in ["partial_3prime", "partial(3)", "3prime_partial"]:
            if tag in h_lower:
                return "3prime_partial"
        for tag in ["partial_5prime", "partial(5)", "5prime_partial"]:
            if tag in h_lower:
                return "5prime_partial"
        return "partial"
    # Heuristic fallback: starts with M = has start codon
    if seq.startswith("M"):
        return "complete"
    return "partial"


COMPL_RANK = {"complete": 2, "3prime_partial": 1, "5prime_partial": 1, "partial": 0}


# --- Parse protein and CDS FASTAs ---

prot_entries = read_fasta(gmst_faa)
cds_entries = read_fasta(gmst_fnn)

if len(prot_entries) != len(cds_entries):
    print(f"WARNING: protein ({len(prot_entries)}) and CDS ({len(cds_entries)}) "
          f"entry counts differ — using protein FASTA only", file=sys.stderr)

# Build CDS lookup by header first-token
cds_by_id = {}
for hdr, seq in cds_entries:
    cds_by_id[hdr.split()[0]] = (hdr, seq)

# Group ORFs by gene
genes = defaultdict(list)
for hdr, seq in prot_entries:
    gene_id = extract_gene_id(hdr, known_genes)
    completeness = infer_completeness(hdr, seq)
    orf_token = hdr.split()[0]
    genes[gene_id].append((len(seq), completeness, orf_token, seq))

# Select best per gene: complete > partial, then longest
n_total = sum(len(v) for v in genes.values())

with open(out_faa, "w") as faa, open(out_fnn, "w") as fnn, open(out_map, "w") as mapf:
    mapf.write("gene_id\torf_id\tcompleteness\tprot_length\n")
    for gene_id in sorted(genes):
        orfs = genes[gene_id]
        orfs.sort(key=lambda x: (COMPL_RANK.get(x[1], 0), x[0]), reverse=True)
        best_len, best_compl, best_token, best_seq = orfs[0]

        faa.write(f">{gene_id}\n{best_seq}\n")
        mapf.write(f"{gene_id}\t{best_token}\t{best_compl}\t{best_len}\n")

        # Write matching CDS if available
        if best_token in cds_by_id:
            _, cds_seq = cds_by_id[best_token]
            fnn.write(f">{gene_id}\n{cds_seq}\n")

print(f"GeneMarkS-T: {len(genes)} best proteins from {n_total} total ORFs "
      f"({n_total - len(genes)} secondary removed)")
