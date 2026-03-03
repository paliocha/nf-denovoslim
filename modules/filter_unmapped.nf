/*
 * FILTER_UNMAPPED — Remove proteins from unmapped transcripts lacking expression.
 *
 * Only runs when genome-guided locus clustering is enabled (--reference_genome).
 * For each protein whose source transcript was "unmapped" in LOCUS_CLUSTER,
 * requires TPM >= min_tpm in >= min_samples Salmon samples to retain it.
 * Gene-assigned and intergenic proteins are always kept.
 */

process FILTER_UNMAPPED {
    tag "${species_label}"

    input:
    path(proteins)
    path(locus_map)
    path(quant_dirs)
    val(species_label)

    output:
    path("${species_label}_filtered.faa"), emit: faa
    path("filter_unmapped_stats.txt"),     emit: stats

    script:
    def min_tpm     = params.unmapped_min_tpm     ?: 1.0
    def min_samples = params.unmapped_min_samples  ?: 2
    """
    filter_unmapped.py \\
        --proteins ${proteins} \\
        --locus-map ${locus_map} \\
        --quant-dirs ${quant_dirs} \\
        --min-tpm ${min_tpm} \\
        --min-samples ${min_samples} \\
        --out ${species_label}_filtered.faa \\
        --stats filter_unmapped_stats.txt
    """
}
