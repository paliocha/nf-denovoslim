/*
 * MMseqs2 nucleotide clustering — remove near-identical transcripts
 */

process MMSEQS2_CLUSTER_NT {
    label 'process_medium'
    tag "${trinity_fasta.baseName}"

    input:
    path(trinity_fasta)

    output:
    path("trinity_clu_nt_rep_seq.fasta"), emit: rep_fasta
    path("trinity_clu_nt_cluster.tsv"),   emit: cluster_tsv

    script:
    def args = task.ext.args ?: ''
    """
    mmseqs easy-cluster \\
        ${trinity_fasta} \\
        trinity_clu_nt \\
        tmp_nt \\
        ${args} \\
        --threads ${task.cpus}
    """
}
