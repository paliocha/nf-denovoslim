/*
 * Salmon index — reusable for both initial and final indexing
 */

process SALMON_INDEX_INITIAL {
    label 'process_medium'
    tag "initial"

    input:
    path(fasta)

    output:
    path("salmon_idx"), emit: index

    script:
    """
    salmon index \\
        -t ${fasta} \\
        -i salmon_idx \\
        --keepDuplicates \\
        -p ${task.cpus}
    """
}

process SALMON_INDEX_FINAL {
    label 'process_medium'
    tag "final"

    input:
    path(fasta)

    output:
    path("st_salmon_idx"), emit: index

    script:
    """
    salmon index \\
        -t ${fasta} \\
        -i st_salmon_idx \\
        -p ${task.cpus}
    """
}
