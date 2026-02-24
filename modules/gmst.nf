/*
 * GMST_PREDICT — GeneMarkS-T self-training ab initio gene prediction
 * Lomsadze et al. (2014) NAR 42(12):e119.  Self-trains on the input
 * transcriptome, then predicts ORFs.  Captures genes with no homology
 * to SwissProt and unusual codon usage that PSAURON may score poorly.
 *
 * Outputs best protein + CDS per gene (longest ORF when multiple).
 */

process GMST_PREDICT {
    tag "${species_label}"

    input:
    path(fasta)
    val(species_label)

    output:
    path("gmst_best.faa"),  emit: faa
    path("gmst_best.fnn"),  emit: fnn
    path("gmst_map.tsv"),   emit: map

    script:
    """
    # Run GeneMarkS-T (self-training mode)
    gmst.pl --faa --fnn --format GFF3 ${fasta}

    # Pick best ORF per gene: completeness → length (PSAURON re-ranks in merge)
    gmst_select_best.py \\
        ${fasta} \\
        ${fasta}.faa \\
        ${fasta}.fnn \\
        gmst_best.faa \\
        gmst_best.fnn \\
        gmst_map.tsv
    """
}
