/*
 * TD2 LongOrfs — extract candidate ORFs from SuperTranscripts
 */

process TD2_LONGORFS {
    label 'process_low'
    tag "${params.species_label}"

    // TD2 writes to ${input}.TD2_dir/ beside the input file
    stageInMode 'copy'

    input:
    path(supertranscripts_fasta)

    output:
    path("${supertranscripts_fasta}.TD2_dir/longest_orfs.pep"), emit: longest_orfs_pep
    path("${supertranscripts_fasta}.TD2_dir/longest_orfs.cds"), emit: longest_orfs_cds
    path("${supertranscripts_fasta}.TD2_dir/longest_orfs.gff3"), emit: longest_orfs_gff3
    path("${supertranscripts_fasta}.TD2_dir"),                   emit: td2_dir

    script:
    def strand_flag = params.td2_strand_specific ? '-S' : ''
    """
    TD2.LongOrfs \\
        -t ${supertranscripts_fasta} \\
        ${strand_flag} \\
        -m ${params.td2_min_orf_length}
    """
}
