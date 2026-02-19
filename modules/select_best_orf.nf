/*
 * SELECT_BEST_ORF — one protein per gene, selected by PSAURON score
 */

process SELECT_BEST_ORF {
    tag "${species_label}"

    input:
    path(psauron_scores)
    path(td2_pep)
    path(td2_gff3)
    val(species_label)

    output:
    path("${species_label}.faa"),      emit: faa
    path("best_orfs.gff3"),            emit: gff3
    path("orf_to_gene_map.tsv"),       emit: map

    script:
    """
    select_best_orf.py \\
        ${psauron_scores} \\
        ${td2_pep} \\
        ${td2_gff3} \\
        ${species_label}
    """
}
