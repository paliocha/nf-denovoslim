/*
 * Lace 2.0 — build one SuperTranscript per gene from Corset clusters
 */

process LACE {
    tag "${species_label}"

    input:
    path(trinity_fasta)
    path(corset_clusters)
    val(species_label)

    output:
    path("supertranscripts.fasta"), emit: fasta

    script:
    """
    export PYTHONUNBUFFERED=1

    Lace \\
        ${trinity_fasta} \\
        ${corset_clusters} \\
        -t \\
        --cores ${task.cpus} \\
        -o lace_out

    cp lace_out/SuperDuper.fasta supertranscripts.fasta
    """
}
