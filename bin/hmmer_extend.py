#!/usr/bin/env python3
"""Pfam domain-guided protein extension and rescue using pyhmmer.

Replaces the merged protein set with an improved version where:
  - Truncated proteins are EXTENDED if Pfam domain evidence supports
    a longer ORF (>= 1.1x current length) on the same transcript.
  - Genes WITHOUT any protein are RESCUED if a Pfam domain hit
    identifies a credible ORF (>= 50 aa) on the transcript.

All processing is in-memory via pyhmmer — no intermediate files,
no text parsing, no hmmpress required.

Usage:
    hmmer_extend.py --merged-faa X --transcripts X --pfam-hmm X \
                    --species X [--cpus 32]

Outputs:
    {species}.extended.faa    — extended + rescued protein FASTA
    hmmer_extend_stats.txt    — summary statistics
    hmmer_extend_map.tsv      — per-gene action log
"""

import argparse
import sys
from collections import defaultdict

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMFile
from pyhmmer.hmmer import hmmsearch


# ── Standard genetic code ─────────────────────────────────────────────────

_CODONS = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

_COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNN')


def _revcomp(seq):
    return seq[::-1].translate(_COMP)


def _translate(nuc):
    """Translate nucleotide string -> protein string (standard code)."""
    prot = []
    for i in range(0, len(nuc) - 2, 3):
        prot.append(_CODONS.get(nuc[i:i + 3].upper(), 'X'))
    return ''.join(prot)


# ── FASTA I/O ─────────────────────────────────────────────────────────────

def parse_fasta(path):
    """Parse FASTA file -> dict {header_id: sequence}."""
    seqs = {}
    cur_id = None
    cur_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if cur_id is not None:
                    seqs[cur_id] = ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id is not None:
        seqs[cur_id] = ''.join(cur_seq)
    return seqs


# ── 6-frame translation ──────────────────────────────────────────────────

def sixframe_orfs(gene_id, nuc_seq, min_len=30):
    """Yield (orf_name, protein_seq) for all ORFs >= min_len aa in 6 frames.

    ORF names encode origin for back-mapping::

        {gene_id}::f{frame}::{aa_offset}

    where aa_offset is the position in the frame's full translation where
    this ORF begins (after splitting at stop codons).
    """
    nuc_seq = nuc_seq.upper().replace('N', 'A')  # mask ambiguous bases
    rc_seq = _revcomp(nuc_seq)

    for frame in range(6):
        strand_seq = nuc_seq if frame < 3 else rc_seq
        offset = frame % 3
        trimmed = strand_seq[offset:]
        trim_len = len(trimmed) - (len(trimmed) % 3)
        if trim_len < 3:
            continue
        trimmed = trimmed[:trim_len]
        protein = _translate(trimmed)

        # Split at stop codons and yield ORFs
        pos = 0
        for orf in protein.split('*'):
            if len(orf) >= min_len:
                yield f"{gene_id}::f{frame}::{pos}", orf
            pos += len(orf) + 1  # +1 for the stop codon


# ── Domain hit processing ────────────────────────────────────────────────

def find_best_orf(gene_hits, orf_seqs, evalue_threshold=1e-5):
    """Select the best ORF for a gene based on Pfam domain evidence.

    Groups domain hits by ORF, merges overlapping domain envelopes to
    compute total covered residues, and picks the ORF with greatest
    coverage.  The returned ORF is trimmed to the nearest upstream Met.

    Returns
    -------
    (protein_seq, pfam_domain_name, domain_coverage_fraction)
    or (None, None, 0.0) if no usable ORF is found.
    """
    good = [h for h in gene_hits if h['i_evalue'] <= evalue_threshold]
    if not good:
        return None, None, 0.0

    # Group by ORF
    by_orf = defaultdict(list)
    for h in good:
        by_orf[h['orf_name']].append(h)

    best_orf_name = None
    best_covered = 0
    best_domain = None

    for orf_name, hits in by_orf.items():
        # Merge overlapping domain envelopes
        intervals = sorted((h['env_from'], h['env_to']) for h in hits)
        merged = []
        for s, e in intervals:
            if merged and s <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        covered = sum(e - s + 1 for s, e in merged)

        if covered > best_covered:
            best_covered = covered
            best_orf_name = orf_name
            best_domain = min(hits, key=lambda h: h['i_evalue'])['domain']

    if best_orf_name is None:
        return None, None, 0.0

    orf_seq = orf_seqs.get(best_orf_name)
    if not orf_seq:
        return None, None, 0.0

    # Trim to nearest upstream Met for a proper start codon
    min_env = min(h['env_from'] for h in by_orf[best_orf_name])
    m_pos = -1
    for i in range(min_env - 1, -1, -1):   # 1-based -> 0-based search
        if orf_seq[i] == 'M':
            m_pos = i
            break

    if m_pos >= 0:
        orf_seq = orf_seq[m_pos:]

    coverage = best_covered / len(orf_seq) if orf_seq else 0.0
    return orf_seq, best_domain, coverage


# ── Helpers ───────────────────────────────────────────────────────────────

def _median(vals):
    if not vals:
        return 0
    s = sorted(vals)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Pfam domain-guided protein extension/rescue (pyhmmer)')
    ap.add_argument('--merged-faa', required=True,
                    help='Merged protein FASTA from merge_predictions')
    ap.add_argument('--transcripts', required=True,
                    help='Corrected representative transcripts FASTA')
    ap.add_argument('--pfam-hmm', required=True,
                    help='Path to Pfam-A.hmm (plain or hmmpress\'d)')
    ap.add_argument('--species', required=True)
    ap.add_argument('--cpus', type=int, default=0,
                    help='CPUs for hmmsearch (0 = all available)')
    ap.add_argument('--min-extend-ratio', type=float, default=1.1,
                    help='Replace protein if domain ORF >= ratio x current '
                         '(default: 1.1)')
    ap.add_argument('--min-rescue-len', type=int, default=50,
                    help='Min ORF length (aa) for rescue (default: 50)')
    ap.add_argument('--evalue', type=float, default=1e-5,
                    help='Domain i-Evalue threshold for extension/rescue '
                         '(default: 1e-5)')
    ap.add_argument('--search-evalue', type=float, default=1e-3,
                    help='E-value pre-filter for hmmsearch (default: 1e-3)')
    args = ap.parse_args()

    # ── Load merged proteins ──────────────────────────────────────────
    merged = parse_fasta(args.merged_faa)
    print(f"Loaded {len(merged):,} merged proteins", file=sys.stderr)

    # ── Load transcripts and 6-frame translate ────────────────────────
    transcripts = parse_fasta(args.transcripts)
    print(f"Loaded {len(transcripts):,} transcripts", file=sys.stderr)

    orf_seqs = {}                      # orf_name -> protein_seq
    gene_orfs = defaultdict(list)      # gene_id -> [orf_names]
    n_orfs = 0

    for gene_id, nuc_seq in transcripts.items():
        for orf_name, orf_seq in sixframe_orfs(gene_id, nuc_seq, min_len=30):
            orf_seqs[orf_name] = orf_seq
            gene_orfs[gene_id].append(orf_name)
            n_orfs += 1

    print(f"6-frame translation: {n_orfs:,} ORFs from "
          f"{len(gene_orfs):,} genes", file=sys.stderr)

    # ── Build pyhmmer digital target database ─────────────────────────
    alphabet = Alphabet.amino()
    digital_targets = []
    for orf_name, orf_seq in orf_seqs.items():
        ts = TextSequence(name=orf_name, sequence=orf_seq)
        digital_targets.append(ts.digitize(alphabet))

    print(f"Prepared {len(digital_targets):,} digital sequences",
          file=sys.stderr)

    # ── hmmsearch: stream Pfam profiles against 6-frame ORFs ──────────
    gene_hits = defaultdict(list)      # gene_id -> [hit_dicts]
    n_hmms = 0

    print("Running hmmsearch against Pfam-A...", file=sys.stderr)

    with HMMFile(args.pfam_hmm) as hmm_file:
        for top_hits in hmmsearch(hmm_file, digital_targets,
                                  cpus=args.cpus,
                                  E=args.search_evalue,
                                  domE=args.search_evalue):
            n_hmms += 1
            if n_hmms % 2000 == 0:
                print(f"  {n_hmms:,} profiles searched...",
                      file=sys.stderr)

            for hit in top_hits:
                if not hit.included:
                    continue
                orf_name = hit.name if isinstance(hit.name, str) else hit.name.decode()
                gene_id = orf_name.split('::')[0]
                qname = top_hits.query_name if isinstance(top_hits.query_name, str) else top_hits.query_name.decode()

                for domain in hit.domains:
                    if not domain.included:
                        continue
                    gene_hits[gene_id].append({
                        'orf_name': orf_name,
                        'domain': qname,
                        'env_from': domain.env_from,
                        'env_to': domain.env_to,
                        'i_evalue': domain.i_evalue,
                    })

    print(f"hmmsearch done: {n_hmms:,} profiles, "
          f"{len(gene_hits):,} genes with domain hits", file=sys.stderr)

    # ── Extension and rescue ──────────────────────────────────────────
    n_extended = 0
    n_rescued = 0
    n_kept = 0
    n_no_protein = 0
    extended_details = []          # (old_len, new_len)
    rescued_lengths = []

    with open(f"{args.species}.extended.faa", 'w') as faa_out, \
         open("hmmer_extend_map.tsv", 'w') as map_out:

        map_out.write("gene_id\taction\told_length\tnew_length\tdomain\n")

        for gene_id in sorted(transcripts.keys()):
            existing = merged.get(gene_id)
            hits = gene_hits.get(gene_id, [])

            if existing and not hits:
                # No domain evidence -> keep as-is
                faa_out.write(f">{gene_id}\n{existing}\n")
                n_kept += 1
                map_out.write(f"{gene_id}\tkept\t{len(existing)}\t"
                              f"{len(existing)}\t-\n")

            elif existing and hits:
                # Try to extend
                orf_seq, domain, _ = find_best_orf(
                    hits, orf_seqs, args.evalue)

                if (orf_seq and
                        len(orf_seq) >= len(existing) * args.min_extend_ratio):
                    faa_out.write(f">{gene_id}\n{orf_seq}\n")
                    n_extended += 1
                    extended_details.append((len(existing), len(orf_seq)))
                    map_out.write(f"{gene_id}\textended\t{len(existing)}\t"
                                  f"{len(orf_seq)}\t{domain}\n")
                else:
                    faa_out.write(f">{gene_id}\n{existing}\n")
                    n_kept += 1
                    map_out.write(f"{gene_id}\tkept\t{len(existing)}\t"
                                  f"{len(existing)}\t-\n")

            elif not existing and hits:
                # Rescue — gene had no protein
                orf_seq, domain, _ = find_best_orf(
                    hits, orf_seqs, args.evalue)

                if orf_seq and len(orf_seq) >= args.min_rescue_len:
                    faa_out.write(f">{gene_id}\n{orf_seq}\n")
                    n_rescued += 1
                    rescued_lengths.append(len(orf_seq))
                    map_out.write(f"{gene_id}\trescued\t0\t"
                                  f"{len(orf_seq)}\t{domain}\n")
                else:
                    n_no_protein += 1

            else:
                # No protein, no hits -> nothing to do
                n_no_protein += 1

    # ── Statistics ────────────────────────────────────────────────────
    total = n_kept + n_extended + n_rescued

    stats = f"""HMMER protein extension statistics (pyhmmer)
{'=' * 50}
Input merged proteins:       {len(merged):>10,}
Total representative genes:  {len(transcripts):>10,}
6-frame ORFs generated:      {n_orfs:>10,}
Pfam profiles searched:      {n_hmms:>10,}
Genes with Pfam domain hits: {len(gene_hits):>10,}

Actions:
  Kept (unchanged):           {n_kept:>10,}
  Extended (domain-guided):   {n_extended:>10,}
  Rescued (new protein):      {n_rescued:>10,}
  No protein (no evidence):   {n_no_protein:>10,}
  Total output proteins:      {total:>10,}
"""

    if extended_details:
        old_lens = [x[0] for x in extended_details]
        new_lens = [x[1] for x in extended_details]
        ratios = [n / o for o, n in extended_details if o > 0]
        stats += f"""
Extension details:
  Old length: mean={sum(old_lens) / len(old_lens):.0f}  median={_median(old_lens)}
  New length: mean={sum(new_lens) / len(new_lens):.0f}  median={_median(new_lens)}
  Ratio: mean={sum(ratios) / len(ratios):.2f}  median={_median(ratios):.2f}
"""

    if rescued_lengths:
        stats += f"""
Rescue details:
  Length: mean={sum(rescued_lengths) / len(rescued_lengths):.0f}  median={_median(rescued_lengths)}
"""

    with open("hmmer_extend_stats.txt", 'w') as f:
        f.write(stats)

    print(stats, file=sys.stderr)


if __name__ == '__main__':
    main()
