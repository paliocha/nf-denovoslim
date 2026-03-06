/*
 * BUILD_TX2GENE — build transcript-to-gene mapping for tximport.
 *
 * Chains Corset clusters → nt dedup → locus clustering → protein dedup
 * to produce a two-column (TXNAME, GENEID) map compatible with
 * tximport::tximport() in R.
 *
 * Quantification against the full decontaminated Trinity assembly
 * plus this tx2gene map gives gene-level counts that include multi-
 * mapping resolution across all isoforms.
 */

process BUILD_TX2GENE {
    tag "${species_label}"

    input:
    path(corset_clusters)
    path(decontam_fasta)
    path(nt_dedup_tsv)
    path(locus_map)
    path(protein_dedup_map)
    val(species_label)

    output:
    path("tx2gene.tsv"),       emit: tx2gene
    path("tx2gene_stats.txt"), emit: stats

    script:
    def locus_arg   = locus_map.name   != 'NO_FILE' ? "--locus-map ${locus_map}" : ''
    def protein_arg = protein_dedup_map.name != 'NO_FILE2' ? "--protein-dedup-map ${protein_dedup_map}" : ''
    """
    build_tx2gene.py \\
        --corset-clusters ${corset_clusters} \\
        --decontam-fasta  ${decontam_fasta} \\
        --nt-dedup-tsv    ${nt_dedup_tsv} \\
        ${locus_arg} \\
        ${protein_arg} \\
        --output tx2gene.tsv \\
        --stats  tx2gene_stats.txt
    """
}
