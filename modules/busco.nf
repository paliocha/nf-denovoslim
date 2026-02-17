/*
 * BUSCO v6 — assess completeness of the final protein set
 */

process BUSCO_QC {
    label 'process_high'
    tag "${species_label}"

    input:
    path(faa)
    val(species_label)

    output:
    path("busco_final"),              emit: outdir
    path("busco_final/short_summary*"), emit: summary

    script:
    def args = task.ext.args ?: ''
    """
    busco \\
        -i ${faa} \\
        ${args} \\
        -m proteins \\
        -o busco_final \\
        -c ${task.cpus}
    """
}
