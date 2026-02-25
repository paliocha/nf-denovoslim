/*
 * DECONTAMINATE_TRINITY — extract original Trinity transcripts from
 * taxonomy-clean Corset clusters.
 *
 * Uses the MMSEQS2_TAXONOMY-filtered representatives to identify which
 * clusters passed the plant filter, then pulls ALL Trinity isoforms
 * from those clusters.  Produces a full decontaminated Trinity FASTA
 * suitable for re-assembly, re-quantification, or sharing.
 */

process DECONTAMINATE_TRINITY {
    tag "${species_label}"

    input:
    path(trinity_fasta)
    path(corset_clusters)
    path(filtered_reps)
    val(species_label)

    output:
    path("${species_label}_trinity_decontaminated.fasta"), emit: fasta
    path("decontamination_stats.txt"),                     emit: stats

    script:
    """
    decontaminate_trinity.py \\
        ${trinity_fasta} \\
        ${corset_clusters} \\
        ${filtered_reps} \\
        ${species_label}_trinity_decontaminated.fasta \\
        decontamination_stats.txt
    """
}
