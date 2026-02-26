/*
 * HMMER_EXTEND — Pfam domain-guided protein extension and rescue (pyhmmer).
 *
 * Single in-memory workflow via pyhmmer (no intermediate files):
 * 1.  6-frame translate all corrected representative transcripts.
 * 2.  hmmsearch Pfam-A.hmm (~20K profiles) against 6-frame ORFs.
 * 3.  Extend truncated proteins where domain evidence supports a
 *     longer ORF (>= 1.1× current length).
 * 4.  Rescue genes without proteins that have Pfam domain hits
 *     with ORF >= 50 aa.
 *
 * Sits between MERGE_PREDICTIONS and MMSEQS2_CLUSTER_PROTEIN.
 * No hmmpress required — pyhmmer reads plain .hmm files.
 */

process HMMER_EXTEND {
    tag "${species_label}"

    input:
    path(merged_faa)
    path(transcripts)
    val(pfam_hmm)
    val(species_label)

    output:
    path("${species_label}.extended.faa"), emit: faa
    path("hmmer_extend_stats.txt"),        emit: stats
    path("hmmer_extend_map.tsv"),          emit: map

    script:
    """
    hmmer_extend.py \\
        --merged-faa ${merged_faa} \\
        --transcripts ${transcripts} \\
        --pfam-hmm ${pfam_hmm} \\
        --species ${species_label} \\
        --cpus ${task.cpus}
    """
}
