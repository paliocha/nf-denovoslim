/*
 * Lace — build SuperTranscripts from Corset clusters
 *
 * Lace constructs a single linear representation (SuperTranscript) per gene
 * by aligning transcripts within each cluster using BLAT and merging them
 * into a non-redundant consensus sequence.
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
    LOCAL=\${TMPDIR:-/tmp}/lace_\$\$
    mkdir -p \$LOCAL
    export MPLCONFIGDIR=\$(mktemp -d -p /tmp)
    export PYTHONUNBUFFERED=1

    Lace \\
        ${deduped_fasta} \\
        ${corset_clusters} \\
        -t \\
        --cores ${task.cpus} \\
        -o \$LOCAL/lace_out

    cp \$LOCAL/lace_out/SuperDuper.fasta supertranscripts.fasta
    rm -rf \$LOCAL
    """
}
