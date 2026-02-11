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
    """
    mmseqs easy-cluster \\
        ${trinity_fasta} \\
        trinity_clu_nt \\
        tmp_nt \\
        --min-seq-id ${params.mmseqs2_nt_id} \\
        -c ${params.mmseqs2_nt_cov} \\
        --cov-mode 1 \\
        --threads ${task.cpus}
    """
}
