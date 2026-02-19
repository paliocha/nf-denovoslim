/*
 * Salmon quantification
 */

process SALMON_QUANT {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)
    path(index)
    val(suffix)

    output:
    path("${sample_id}_${suffix}"), emit: quant_dir

    script:
    def args = task.ext.args ?: ''
    """
    salmon quant \\
        -i ${index} \\
        -l A \\
        -1 ${reads_1} -2 ${reads_2} \\
        --validateMappings \\
        --seqBias --gcBias --posBias \\
        ${args} \\
        -p ${task.cpus} \\
        -o ${sample_id}_${suffix}
    """
}
