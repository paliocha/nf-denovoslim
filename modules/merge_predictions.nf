/*
 * MERGE_PREDICTIONS — merge TD2 + MetaEuk + GeneMarkS-T into a single
 * best-per-gene protein set.  Ranks: completeness → length → PSAURON.
 */

process MERGE_PREDICTIONS {
    tag "${species_label}"

    input:
    path(td2_faa)
    path(td2_map)
    path(metaeuk_faa)
    path(metaeuk_map)
    path(metaeuk_psauron)
    path(gmst_faa)
    path(gmst_map)
    path(gmst_psauron)
    path(salmon_quant_dirs)
    val(min_psauron)
    val(species_label)

    output:
    path("${species_label}.faa"), emit: faa
    path("merge_map.tsv"),        emit: map
    path("merge_stats.txt"),      emit: stats

    script:
    def quant_dirs = salmon_quant_dirs.collect { d -> d.name }.join(',')
    """
    merge_predictions.py \\
        --td2-faa ${td2_faa} \\
        --td2-map ${td2_map} \\
        --metaeuk-faa ${metaeuk_faa} \\
        --metaeuk-map ${metaeuk_map} \\
        --metaeuk-psauron ${metaeuk_psauron} \\
        --gmst-faa ${gmst_faa} \\
        --gmst-map ${gmst_map} \\
        --gmst-psauron ${gmst_psauron} \\
        --salmon-quants ${quant_dirs} \\
        --min-psauron ${min_psauron} \\
        --species ${species_label}
    """
}
