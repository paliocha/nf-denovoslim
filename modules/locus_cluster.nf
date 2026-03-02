/*
 * LOCUS_CLUSTER — collapse transcripts that map to the same genomic locus.
 *
 * Reads minimap2 PAF alignments, groups transcripts by overlapping genomic
 * coordinates (same chromosome, same strand), picks the best representative
 * per locus (most aligned bases), and retains all unmapped transcripts.
 *
 * Sits between MMSEQS2_CLUSTER (nt dedup) and DIAMOND_BLASTX (frameshift
 * correction).  Only runs when --reference_genome is provided.
 */

process LOCUS_CLUSTER {
    tag "${species_label}"

    input:
    path(query_fasta)
    path(paf)
    val(species_label)

    output:
    path("${species_label}_locus_collapsed.fasta"), emit: fasta
    path("locus_map.tsv"),                          emit: map
    path("locus_stats.txt"),                        emit: stats

    script:
    def max_intron   = params.locus_max_intron   ?: 200000
    def min_coverage = params.locus_min_coverage ?: 0.5
    def min_identity = params.locus_min_identity ?: 0.7
    """
    locus_cluster.py \\
        --paf ${paf} \\
        --fasta ${query_fasta} \\
        --out ${species_label}_locus_collapsed.fasta \\
        --map locus_map.tsv \\
        --stats locus_stats.txt \\
        --max-intron ${max_intron} \\
        --min-coverage ${min_coverage} \\
        --min-identity ${min_identity}
    """
}
