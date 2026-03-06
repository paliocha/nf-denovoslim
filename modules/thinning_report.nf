/*
 * THINNING_REPORT — pipeline summary statistics
 */

process THINNING_REPORT {
    tag "${species_label}"

    input:
    path(trinity_fasta)
    path(representatives_fasta)
    path(corset_clust)
    path(merge_map)
    path(faa)
    path(initial_quant_dirs)
    path(final_quant_dirs)
    path(full_quant_dirs)
    path(tx2gene)
    path(busco_trinity_summary)
    path(busco_reps_summary)
    path(busco_corrected_summary)
    path(busco_final_summary)
    path(id_validation)
    path(sortmerna_logs, stageAs: 'sortmerna_??.log')
    path(taxonomy_breakdown)
    path(transannot_tsv)
    path(merge_stats)
    path(dedup_stats)
    path(protein_dedup_stats)
    val(species_label)

    output:
    path("${species_label}_thinning_report.txt"), emit: report

    script:
    def init_dirs  = initial_quant_dirs.collect { d -> d.name }.join(',')
    def final_dirs = final_quant_dirs.collect { d -> d.name }.join(',')
    def full_dirs  = full_quant_dirs.collect { d -> d.name }.join(',')
    def log_files  = sortmerna_logs.collect { f -> f.name }.join(',')
    """
    thinning_report.py \\
        --species          ${species_label} \\
        --trinity          ${trinity_fasta} \\
        --representatives  ${representatives_fasta} \\
        --clusters         ${corset_clust} \\
        --merge-map        ${merge_map} \\
        --faa              ${faa} \\
        --initial-quants   ${init_dirs} \\
        --final-quants     ${final_dirs} \\
        --full-quants      ${full_dirs} \\
        --tx2gene          ${tx2gene} \\
        --busco-trinity    ${busco_trinity_summary} \\
        --busco-reps       ${busco_reps_summary} \\
        --busco-corrected  ${busco_corrected_summary} \\
        --busco-final      ${busco_final_summary} \\
        --id-validation    ${id_validation} \\
        --sortmerna-logs   ${log_files} \\
        --taxonomy         ${taxonomy_breakdown} \\
        --transannot       ${transannot_tsv} \\
        --merge-stats      ${merge_stats} \\
        --dedup-stats      ${dedup_stats} \\
        --protein-dedup-stats ${protein_dedup_stats}
    """
}
