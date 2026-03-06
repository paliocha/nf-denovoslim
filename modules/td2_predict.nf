/*
 * TD2 Predict — ORF prediction with homology support
 */

process TD2_PREDICT {
    tag "${species_label}"

    // TD2 writes outputs beside the input fasta
    stageInMode 'copy'

    input:
    path(representatives_fasta)
    path(swissprot_m8)
    path(pfam_m8)
    path(td2_dir)
    val(species_label)

    output:
    path("${representatives_fasta}.TD2.pep"),              emit: pep
    path("${representatives_fasta}.TD2.cds"),              emit: cds
    path("${representatives_fasta}.TD2.gff3"),             emit: gff3
    path("${representatives_fasta}.TD2.bed"),              emit: bed
    path("${representatives_fasta.baseName}/psauron_score.csv"), emit: psauron_scores

    script:
    def args = task.ext.args ?: ''
    """
    # Concatenate homology hits from both searches
    cat ${swissprot_m8} ${pfam_m8} > combined_alnRes.m8

    # Ensure TD2 output dir (stem) is present beside the input fasta.
    # Nextflow stages td2_dir with its original name; rename if needed.
    if [ ! -d "${representatives_fasta.baseName}" ] && [ -d "${td2_dir}" ]; then
        mv ${td2_dir} ${representatives_fasta.baseName}
    fi

    TD2.Predict \\
        -t ${representatives_fasta} \\
        --retain-mmseqs-hits combined_alnRes.m8 \\
        ${args} \\
        -v
    """
}
