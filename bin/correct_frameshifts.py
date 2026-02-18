#!/usr/bin/env python3
"""Correct frameshifts in nucleotide sequences using Diamond blastx BTOP strings.

Diamond blastx with -F 15 performs frameshift-tolerant protein alignment.
The BTOP string records frameshift events:
    \\  = +1 shift (extra base inserted in query → delete 1 base to fix)
    /  = -1 shift (base missing from query → insert 1 N to fix)

Usage:
    correct_frameshifts.py <input.fasta> <diamond_fs.tsv> <output.fasta>

The Diamond TSV must have columns (outfmt 6):
    qseqid qstart qend qlen qframe btop
"""

import re
import sys

# --------------------------------------------------------------------------- #
#  BTOP tokenizer
# --------------------------------------------------------------------------- #

# Tokens in a BTOP string:
#   number       → N matching positions
#   XY           → mismatch (query aa X, subject aa Y)
#   -X           → gap in query (insertion in subject)
#   X-           → gap in subject (insertion in query)
#   \            → +1 frameshift
#   /            → -1 frameshift
BTOP_RE = re.compile(r'(\d+|[A-Z*]-|-[A-Z*]|[A-Z*]{2}|[/\\])')


def tokenize_btop(btop_str):
    """Split a BTOP string into individual tokens."""
    return BTOP_RE.findall(btop_str)


# --------------------------------------------------------------------------- #
#  Parse Diamond output → frameshift positions
# --------------------------------------------------------------------------- #

def parse_frameshifts(diamond_tsv):
    """Parse Diamond blastx output and return frameshift positions per query.

    Returns dict: qseqid → list of (position_0based, 'delete'|'insert')
    sorted by position ascending.

    Only considers the best hit (first line) per query on + strand.
    """
    best_hits = {}  # qseqid → (bitscore_proxy, line_fields)

    with open(diamond_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            qseqid = fields[0]
            qframe = int(fields[4])

            # Only consider forward-strand hits (frames 1, 2, 3)
            # SuperTranscripts from Lace should be in the correct orientation
            if qframe < 0:
                continue

            # Keep first hit per query (Diamond outputs best first by default)
            if qseqid not in best_hits:
                best_hits[qseqid] = fields

    frameshifts = {}

    for qseqid, fields in best_hits.items():
        qstart = int(fields[1])  # 1-based
        qend   = int(fields[2])  # 1-based
        qframe = int(fields[4])
        btop   = fields[5]

        tokens = tokenize_btop(btop)

        # Check if there are any frameshifts at all
        if '\\' not in tokens and '/' not in tokens:
            continue

        # Walk through BTOP, tracking position in query (0-based nucleotide)
        # qstart is 1-based, convert to 0-based
        qpos = qstart - 1
        corrections = []

        for token in tokens:
            if token.isdigit():
                # N matching amino acid positions → 3*N query bases consumed
                qpos += 3 * int(token)
            elif token == '\\':
                # +1 frameshift: 1 extra base at qpos → delete it
                corrections.append((qpos, 'delete'))
                qpos += 1
            elif token == '/':
                # -1 frameshift: base missing at qpos → insert N
                corrections.append((qpos, 'insert'))
                qpos -= 1
            elif token.startswith('-'):
                # Gap in query (deletion in DNA vs protein): 0 bases consumed
                pass
            elif token.endswith('-'):
                # Gap in subject (insertion in DNA vs protein): 3 bases consumed
                qpos += 3
            elif len(token) == 2 and token[0].isalpha() and token[1].isalpha():
                # Mismatch: 3 bases consumed
                qpos += 3
            # else: unknown token, skip

        if corrections:
            # Sort by position (should already be in order, but ensure it)
            corrections.sort(key=lambda x: x[0])
            frameshifts[qseqid] = corrections

    return frameshifts


# --------------------------------------------------------------------------- #
#  Apply corrections to sequences
# --------------------------------------------------------------------------- #

def apply_corrections(seq, corrections):
    """Apply frameshift corrections to a nucleotide sequence.

    Corrections are applied from END to START to preserve upstream positions.
    """
    seq_list = list(seq)

    # Work backwards through corrections
    for pos, action in reversed(corrections):
        if pos < 0 or pos > len(seq_list):
            continue
        if action == 'delete' and pos < len(seq_list):
            del seq_list[pos]
        elif action == 'insert':
            seq_list.insert(pos, 'N')

    return ''.join(seq_list)


# --------------------------------------------------------------------------- #
#  Main: read FASTA, apply corrections, write output
# --------------------------------------------------------------------------- #

def read_fasta(path):
    """Simple FASTA reader returning list of (header, sequence) tuples."""
    sequences = []
    header = None
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    sequences.append((header, ''.join(seq_parts)))
                header = line[1:].split()[0]  # first word only
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            sequences.append((header, ''.join(seq_parts)))
    return sequences


def write_fasta(path, sequences, line_width=80):
    """Write sequences as FASTA with fixed line width."""
    with open(path, 'w') as fh:
        for header, seq in sequences:
            fh.write(f'>{header}\n')
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i+line_width] + '\n')


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.fasta> <diamond_fs.tsv> <output.fasta>",
              file=sys.stderr)
        sys.exit(1)

    input_fasta  = sys.argv[1]
    diamond_tsv  = sys.argv[2]
    output_fasta = sys.argv[3]

    # Parse frameshifts from Diamond output
    frameshifts = parse_frameshifts(diamond_tsv)
    n_with_hits = len(frameshifts)  # sequences that had frameshift-containing hits

    # Read input sequences
    sequences = read_fasta(input_fasta)

    # Apply corrections
    n_corrected = 0
    n_total_fs = 0
    n_deletions = 0
    n_insertions = 0
    corrected_sequences = []

    for header, seq in sequences:
        if header in frameshifts:
            corrections = frameshifts[header]
            corrected_seq = apply_corrections(seq, corrections)
            corrected_sequences.append((header, corrected_seq))
            n_corrected += 1
            n_total_fs += len(corrections)
            n_deletions += sum(1 for _, a in corrections if a == 'delete')
            n_insertions += sum(1 for _, a in corrections if a == 'insert')
        else:
            corrected_sequences.append((header, seq))

    # Write output
    write_fasta(output_fasta, corrected_sequences)

    # Print stats to stdout
    print(f"Frameshift correction summary:")
    print(f"  Total sequences:          {len(sequences)}")
    print(f"  Sequences with frameshifts: {n_with_hits}")
    print(f"  Sequences corrected:      {n_corrected}")
    print(f"  Total frameshifts fixed:  {n_total_fs}")
    print(f"    Deletions (extra base): {n_deletions}")
    print(f"    Insertions (N added):   {n_insertions}")
    print(f"  Sequences unchanged:      {len(sequences) - n_corrected}")


if __name__ == '__main__':
    main()
