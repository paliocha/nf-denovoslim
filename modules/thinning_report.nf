/*
 * THINNING_REPORT — pipeline summary statistics
 */

process THINNING_REPORT {
    tag "${species_label}"

    input:
    path(trinity_fasta)
    path(supertranscripts_fasta)
    path(corset_clust)
    path(orf_to_gene_map)
    path(faa)
    path(initial_quant_dirs)
    path(final_quant_dirs)
    path(busco_trinity_summary)
    path(busco_final_summary)
    path(id_validation)
    path(sortmerna_logs, stageAs: 'sortmerna_??.log')
    path(taxonomy_breakdown)
    val(species_label)

    output:
    path("${species_label}_thinning_report.txt"), emit: report

    script:
    def init_dirs  = initial_quant_dirs.collect { d -> d.name }.join(',')
    def final_dirs = final_quant_dirs.collect { d -> d.name }.join(',')
    def log_files  = sortmerna_logs.collect { f -> f.name }.join(',')
    """
    thinning_report.py \\
        ${species_label} \\
        ${trinity_fasta} \\
        ${supertranscripts_fasta} \\
        ${corset_clust} \\
        ${orf_to_gene_map} \\
        ${faa} \\
        ${init_dirs} \\
        ${final_dirs} \\
        ${busco_trinity_summary} \\
        ${busco_final_summary} \\
        ${id_validation} \\
        ${log_files} \\
        ${taxonomy_breakdown}
    """
}
