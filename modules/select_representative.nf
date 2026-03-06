/*
 * SELECT_REP — pick the longest Trinity transcript per Corset cluster
 */

process SELECT_REP {
    tag "${species_label}"

    input:
    path(trinity_fasta)
    path(corset_clusters)
    val(species_label)

    output:
    path("representatives.fasta"), emit: fasta
    path("rep_map.tsv"),           emit: map

    script:
    """
    select_representative.py \\
        ${trinity_fasta} \\
        ${corset_clusters} \\
        representatives.fasta \\
        rep_map.tsv
    """
}
