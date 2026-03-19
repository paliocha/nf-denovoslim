/*
 * FILTER_UNMAPPED — Remove proteins from unmapped transcripts lacking evidence.
 *
 * Only runs when genome-guided locus clustering is enabled (--reference_genome).
 * For each protein whose source transcript was "unmapped" in LOCUS_CLUSTER:
 *   1. If it also has no taxonomy hit (taxid=0): always removed (contaminant).
 *   2. Otherwise: requires TPM >= min_tpm in >= min_samples to retain.
 * Gene-assigned and intergenic proteins are always kept.
 */

process FILTER_UNMAPPED {
    tag "${species_label}"

    input:
    path(proteins)
    path(locus_map)
    path(taxonomy_tsv)
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
        --taxonomy-tsv ${taxonomy_tsv} \\
        --quant-dirs ${quant_dirs} \\
        --min-tpm ${min_tpm} \\
        --min-samples ${min_samples} \\
        --out ${species_label}_filtered.faa \\
        --stats filter_unmapped_stats.txt
    """
}
