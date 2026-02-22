/*
 * BUSCO v6 — transcriptome or protein completeness assessment
 */

process BUSCO {
    tag "${species_label}_${suffix}"

    input:
    path(fasta)
    val(species_label)
    val(suffix)

    output:
    path("busco_${suffix}"),                   emit: outdir
    path("busco_${suffix}/short_summary*.txt"), emit: summary
    path("busco_${suffix}/short_summary*"),     emit: summary_all

    script:
    def args = task.ext.args ?: ''
    """
    busco \\
        -i ${fasta} \\
        ${args} \\
        -o busco_${suffix} \\
        -c ${task.cpus}
    """
}
