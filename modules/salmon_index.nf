/*
 * Salmon index — builds a Salmon index from a FASTA reference
 *
 * Usage with DSL2 aliasing:
 *   include { SALMON_INDEX as SALMON_INDEX_INITIAL } from './modules/salmon_index'
 *   include { SALMON_INDEX as SALMON_INDEX_FINAL   } from './modules/salmon_index'
 *
 * Per-alias flags via conf/base.config:
 *   withName: 'SALMON_INDEX_INITIAL' { ext.args = '--keepDuplicates' }
 */

process SALMON_INDEX {
    tag "${fasta.simpleName}"

    input:
    path(fasta)

    output:
    path("salmon_idx"), emit: index

    script:
    def args = task.ext.args ?: ''
    """
    salmon index \\
        -t ${fasta} \\
        -i salmon_idx \\
        ${args} \\
        -p ${task.cpus}
    """
}
