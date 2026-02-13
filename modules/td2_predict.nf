/*
 * TD2 Predict — final ORF prediction with homology support
 */

process TD2_PREDICT {
    label 'process_medium'
    tag "${params.species_label}"

    // TD2 writes outputs beside the input fasta
    stageInMode 'copy'

    input:
    path(supertranscripts_fasta)
    path(swissprot_m8)
    path(pfam_m8)
    path(td2_dir)

    output:
    path("${supertranscripts_fasta}.TD2.pep"),              emit: pep
    path("${supertranscripts_fasta}.TD2.cds"),              emit: cds
    path("${supertranscripts_fasta}.TD2.gff3"),             emit: gff3
    path("${supertranscripts_fasta}.TD2.bed"),              emit: bed
    path("${supertranscripts_fasta}.TD2_dir/psauron_score.csv"), emit: psauron_scores

    script:
    """
    # Concatenate homology hits from both searches
    cat ${swissprot_m8} ${pfam_m8} > combined_alnRes.m8

    # Ensure TD2_dir is where TD2.Predict expects it (beside the input fasta)
    ln -sf ${td2_dir} ${supertranscripts_fasta}.TD2_dir

    TD2.Predict \\
        -t ${supertranscripts_fasta} \\
        --retain-mmseqs-hits combined_alnRes.m8
    """
}
