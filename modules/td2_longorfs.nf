/*
 * TD2 LongOrfs — extract candidate ORFs from SuperTranscripts
 */

process TD2_LONGORFS {
    tag "${species_label}"

    // TD2 writes to ./{stem}/ beside the input file (strips .fasta extension)
    stageInMode 'copy'

    input:
    path(supertranscripts_fasta)
    val(species_label)

    output:
    path("${supertranscripts_fasta.baseName}/longest_orfs.pep"), emit: longest_orfs_pep
    path("${supertranscripts_fasta.baseName}/longest_orfs.cds"), emit: longest_orfs_cds
    path("${supertranscripts_fasta.baseName}/longest_orfs.gff3"), emit: longest_orfs_gff3
    path("${supertranscripts_fasta.baseName}"),                   emit: td2_dir

    script:
    def args = task.ext.args ?: ''
    """
    TD2.LongOrfs \\
        -t ${supertranscripts_fasta} \\
        ${args} \\
        --threads ${task.cpus}
    """
}
