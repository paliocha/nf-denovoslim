/*
 * Salmon index
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
