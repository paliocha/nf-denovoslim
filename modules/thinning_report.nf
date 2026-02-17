/*
 * THINNING_REPORT — final summary of assembly reduction statistics
 */

process THINNING_REPORT {
    label 'process_low'
    tag "${species_label}"

    input:
    path(trinity_fasta)
    path(deduped_fasta)
    path(supertranscripts_fasta)
    path(corset_clust)
    path(orf_to_gene_map)
    path(faa)
    path(initial_quant_dirs)
    path(final_quant_dirs)
    path(busco_summary)
    path(id_validation)
    path(sortmerna_logs)
    val(species_label)

    output:
    path("${species_label}_thinning_report.txt"), emit: report

    script:
    def init_dirs  = initial_quant_dirs.collect { it.name }.join(',')
    def final_dirs = final_quant_dirs.collect { it.name }.join(',')
    def log_files  = sortmerna_logs.collect { it.name }.join(',')
    """
    # Locate BUSCO short summary
    BUSCO_FILE=\$(find ${busco_summary} -name 'short_summary*' -type f | head -1)
    if [ -z "\$BUSCO_FILE" ]; then BUSCO_FILE="/dev/null"; fi

    thinning_report.py \\
        ${species_label} \\
        ${trinity_fasta} \\
        ${deduped_fasta} \\
        ${supertranscripts_fasta} \\
        ${corset_clust} \\
        ${orf_to_gene_map} \\
        ${faa} \\
        ${init_dirs} \\
        ${final_dirs} \\
        \$BUSCO_FILE \\
        ${id_validation} \\
        ${log_files}
    """
}
