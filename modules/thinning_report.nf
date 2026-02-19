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
    path(busco_trinity)
    path(busco_final)
    path(id_validation)
    path(sortmerna_logs, stageAs: 'sortmerna_??.log')
    val(species_label)

    output:
    path("${species_label}_thinning_report.txt"), emit: report

    script:
    def init_dirs  = initial_quant_dirs.collect { it.name }.join(',')
    def final_dirs = final_quant_dirs.collect { it.name }.join(',')
    def log_files  = sortmerna_logs.collect { it.name }.join(',')
    """
    # Locate BUSCO short summaries
    BUSCO_TRINITY_FILE=\$(find ${busco_trinity} -name 'short_summary*' -type f | head -1)
    if [ -z "\$BUSCO_TRINITY_FILE" ]; then BUSCO_TRINITY_FILE="/dev/null"; fi

    BUSCO_FINAL_FILE=\$(find ${busco_final} -name 'short_summary*' -type f | head -1)
    if [ -z "\$BUSCO_FINAL_FILE" ]; then BUSCO_FINAL_FILE="/dev/null"; fi

    thinning_report.py \\
        ${species_label} \\
        ${trinity_fasta} \\
        ${supertranscripts_fasta} \\
        ${corset_clust} \\
        ${orf_to_gene_map} \\
        ${faa} \\
        ${init_dirs} \\
        ${final_dirs} \\
        \$BUSCO_TRINITY_FILE \\
        \$BUSCO_FINAL_FILE \\
        ${id_validation} \\
        ${log_files}
    """
}
