/*
 * BUSCO v6 — assess completeness of the final protein set
 */

process BUSCO_QC {
    label 'process_high'
    tag "${params.species_label}"

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path(faa)

    output:
    path("busco_final"),              emit: outdir
    path("busco_final/short_summary*"), emit: summary

    script:
    """
    busco \\
        -i ${faa} \\
        -l ${params.busco_lineage} \\
        -m proteins \\
        -o busco_final \\
        -c ${task.cpus}
    """
}
