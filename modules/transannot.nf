/*
 * TransAnnot — functional annotation (SwissProt + Pfam + eggNOG)
 *
 * Uses paliocha/transannot (update-eggnog-v7 branch) which natively
 * downloads eggNOG 7 protein families and normalizes pipe/hyphen OG names.
 * No curl wrapper or post-hoc AWK fix needed.
 */

process TRANSANNOT {
    tag "${species_label}"

    input:
    path(faa)
    val(species_label)
    val(pfam_db)
    val(eggnog_db)
    val(swissprot_db)

    output:
    path("${species_label}_transannot.tsv"), emit: annotation

    script:
    """
    transannot createquerydb ${faa} queryDB tmp_createdb

    transannot annotate \\
        queryDB \\
        ${pfam_db} \\
        ${eggnog_db} \\
        ${swissprot_db} \\
        ${species_label}_transannot.tsv \\
        tmp_annotate \\
        --no-run-clust \\
        --threads ${task.cpus}
    """
}
