/*
 * LOCUS_CLUSTER — collapse transcripts that map to the same genomic gene/locus.
 *
 * Two modes:
 *   1. Gene-level (--reference_gff provided): assign each mapped transcript
 *      to the best-overlapping reference gene, keep one per gene.
 *   2. Coordinate-overlap fallback (no GFF): merge overlapping alignment
 *      intervals into ad-hoc loci.
 *
 * Unmapped transcripts are always retained.
 * Sits between MMSEQS2_CLUSTER (nt dedup) and DIAMOND_BLASTX (frameshift
 * correction).  Only runs when --reference_genome is provided.
 */

process LOCUS_CLUSTER {
    tag "${species_label}"

    input:
    path(query_fasta)
    path(paf)
    val(reference_gff)
    val(species_label)

    output:
    path("${species_label}_locus_collapsed.fasta"), emit: fasta
    path("locus_map.tsv"),                          emit: map
    path("locus_stats.txt"),                        emit: stats

    script:
    def max_intron   = params.locus_max_intron   ?: 200000
    def min_coverage = params.locus_min_coverage ?: 0.5
    def min_mapq     = params.locus_min_mapq     ?: 5
    def gene_flank   = params.locus_gene_flank   ?: 5000
    def gff_flag     = reference_gff ? "--gff ${reference_gff}" : ''
    """
    locus_cluster.py \\
        --paf ${paf} \\
        --fasta ${query_fasta} \\
        --out ${species_label}_locus_collapsed.fasta \\
        --map locus_map.tsv \\
        --stats locus_stats.txt \\
        --max-intron ${max_intron} \\
        --min-coverage ${min_coverage} \\
        --min-mapq ${min_mapq} \\
        --gene-flank ${gene_flank} \\
        ${gff_flag}
    """
}
