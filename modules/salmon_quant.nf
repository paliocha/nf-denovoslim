/*
 * Salmon quantification — initial (with --dumpEq for Grouper) and final (gene-level)
 */

process SALMON_QUANT_INITIAL {
    label 'process_high'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)
    path(index)

    output:
    path("${sample_id}_quant"), emit: quant_dir

    script:
    """
    salmon quant \\
        -i ${index} \\
        -l A \\
        -1 ${reads_1} -2 ${reads_2} \\
        --validateMappings \\
        --seqBias --gcBias --posBias \\
        --rangeFactorizationBins 4 \\
        --dumpEq --writeOrphanLinks \\
        -p ${task.cpus} \\
        -o ${sample_id}_quant
    """
}

process SALMON_QUANT_FINAL {
    label 'process_high'
    tag "${sample_id}"

    publishDir "${params.outdir}/salmon_final", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)
    path(index)

    output:
    path("${sample_id}_st_quant"), emit: quant_dir

    script:
    """
    salmon quant \\
        -i ${index} \\
        -l A \\
        -1 ${reads_1} -2 ${reads_2} \\
        --validateMappings \\
        --seqBias --gcBias --posBias \\
        --rangeFactorizationBins 4 \\
        -p ${task.cpus} \\
        -o ${sample_id}_st_quant
    """
}
