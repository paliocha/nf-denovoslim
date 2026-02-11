/*
 * Lace — build SuperTranscripts from Corset clusters
 *
 * Lace constructs a single linear representation (SuperTranscript) per gene
 * by aligning transcripts within each cluster using BLAT and merging them
 * into a non-redundant consensus sequence.
 */

process LACE {
    label 'process_high'
    tag "${params.species_label}"

    publishDir "${params.outdir}/supertranscripts", mode: 'copy'

    input:
    path(deduped_fasta)
    path(corset_clusters)

    output:
    path("supertranscripts.fasta"), emit: fasta

    script:
    """
    Lace \\
        ${deduped_fasta} \\
        ${corset_clusters} \\
        -t \\
        --cores ${task.cpus} \\
        -o lace_out

    cp lace_out/SuperDuper.fasta supertranscripts.fasta
    """
}
