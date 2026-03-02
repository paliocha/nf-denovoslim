/*
 * MINIMAP2_SPLICE — splice-aware mapping of transcripts to a reference genome.
 *
 * Uses minimap2 -x splice for transcript-to-genome alignment, producing
 * PAF output for downstream locus clustering.  Suitable for same-species
 * or cross-genus mapping (unmapped transcripts are retained by LOCUS_CLUSTER).
 *
 * The reference genome is passed as val (not path) to avoid staging a
 * multi-GB file into the work directory — container bind mounts make it
 * accessible directly.
 */

process MINIMAP2_SPLICE {
    tag "${species_label}"

    input:
    path(query_fasta)
    val(reference_genome)
    val(species_label)

    output:
    path("alignments.paf"), emit: paf

    script:
    """
    minimap2 \\
        -x splice \\
        --secondary=no \\
        -t ${task.cpus} \\
        "${reference_genome}" \\
        ${query_fasta} \\
        > alignments.paf

    N_ALN=\$(wc -l < alignments.paf)
    N_QUERY=\$(grep -c '^>' ${query_fasta})
    N_MAPPED=\$(awk '{print \$1}' alignments.paf | sort -u | wc -l)
    echo "minimap2 splice: \$N_QUERY queries, \$N_ALN alignments, \$N_MAPPED mapped"
    """
}
