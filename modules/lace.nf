/*
 * Lace — build SuperTranscripts from Corset clusters
 *
 * Lace constructs a single linear representation (SuperTranscript) per gene
 * by aligning transcripts within each cluster using BLAT and merging them
 * into a non-redundant consensus sequence.
 *
 * Lace creates one FASTA file per cluster in its output dir — thousands of
 * tiny files that are very slow on NFS. With scratch enabled (Orion), the
 * task CWD is already on node-local SSD, so this I/O stays local.
 */

process LACE {
    label 'process_high'
    tag "${species_label}"

    input:
    path(deduped_fasta)
    path(corset_clusters)
    val(species_label)

    output:
    path("supertranscripts.fasta"), emit: fasta

    script:
    """
    export MPLCONFIGDIR=\$(mktemp -d)
    export PYTHONUNBUFFERED=1

    Lace \\
        ${deduped_fasta} \\
        ${corset_clusters} \\
        -t \\
        --cores ${task.cpus} \\
        -o lace_out

    cp lace_out/SuperDuper.fasta supertranscripts.fasta
    """
}
