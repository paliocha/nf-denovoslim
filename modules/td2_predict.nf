/*
 * TD2 Predict — ORF prediction with homology support
 */

process TD2_PREDICT {
    tag "${species_label}"

    // TD2 writes outputs beside the input fasta
    stageInMode 'copy'

    input:
    path(supertranscripts_fasta)
    path(swissprot_m8)
    path(pfam_m8)
    path(td2_dir)
    val(species_label)

    output:
    path("${supertranscripts_fasta}.TD2.pep"),              emit: pep
    path("${supertranscripts_fasta}.TD2.cds"),              emit: cds
    path("${supertranscripts_fasta}.TD2.gff3"),             emit: gff3
    path("${supertranscripts_fasta}.TD2.bed"),              emit: bed
    path("${supertranscripts_fasta.baseName}/psauron_score.csv"), emit: psauron_scores

    script:
    """
    # Concatenate homology hits from both searches
    cat ${swissprot_m8} ${pfam_m8} > combined_alnRes.m8

    # Ensure TD2 output dir (stem) is present beside the input fasta.
    # Nextflow stages td2_dir with its original name; rename if needed.
    if [ ! -d "${supertranscripts_fasta.baseName}" ] && [ -d "${td2_dir}" ]; then
        mv ${td2_dir} ${supertranscripts_fasta.baseName}
    fi

    TD2.Predict \\
        -t ${supertranscripts_fasta} \\
        --retain-mmseqs-hits combined_alnRes.m8 \\
        -v
    """
}
