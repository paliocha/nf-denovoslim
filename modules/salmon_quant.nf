/*
 * Salmon quantification — maps reads and quantifies expression
 *
 * Usage with DSL2 aliasing:
 *   include { SALMON_QUANT as SALMON_QUANT_INITIAL } from './modules/salmon_quant'
 *   include { SALMON_QUANT as SALMON_QUANT_FINAL   } from './modules/salmon_quant'
 *
 * Per-alias flags via conf/base.config:
 *   withName: 'SALMON_QUANT_INITIAL' { ext.args = '--hardFilter --dumpEq' }
 *   withName: 'SALMON_QUANT_FINAL'   { ext.args = '--rangeFactorizationBins 4'
 *                                       publishDir = [...] }
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
