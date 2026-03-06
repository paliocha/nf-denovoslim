/*
 * MMSEQS2_CLUSTER — sequence-level dedup at 90% nucleotide identity
 * Catches near-identical transcripts that Corset missed (different read sets,
 * assembly duplicates).  Uses mmseqs easy-cluster on representative sequences.
 */

process MMSEQS2_CLUSTER {
    tag "${species_label}"

    input:
    path(representatives_fasta)
    val(species_label)

    output:
    path("representatives_dedup.fasta"), emit: fasta
    path("nt_cluster_map.tsv"),           emit: cluster_tsv
    path("cluster_stats.txt"),           emit: stats

    script:
    """
    # Run MMseqs2 easy-cluster at 90% nucleotide identity
    mmseqs easy-cluster \\
        ${representatives_fasta} \\
        clust \\
        tmp_clust \\
        --min-seq-id 0.90 \\
        -c 0.8 \\
        --cov-mode 0 \\
        --threads ${task.cpus}

    # Extract representative sequences (cluster representatives)
    # easy-cluster outputs: clust_rep_seq.fasta, clust_cluster.tsv, clust_all_seqs.fasta
    cp clust_rep_seq.fasta representatives_dedup.fasta
    cp clust_cluster.tsv nt_cluster_map.tsv

    # Report
    N_IN=\$(grep -c '^>' ${representatives_fasta})
    N_OUT=\$(grep -c '^>' representatives_dedup.fasta)
    N_COLLAPSED=\$((N_IN - N_OUT))
    PCT=\$(awk "BEGIN{printf \\"%.1f\\", \$N_COLLAPSED/\$N_IN*100}")

    cat > cluster_stats.txt <<EOF
    MMseqs2 nucleotide dedup (90% identity, 80% coverage)
    Input sequences:     \$N_IN
    Output sequences:    \$N_OUT
    Collapsed:           \$N_COLLAPSED (\${PCT}%)
    EOF

    echo "MMseqs2 dedup: \$N_IN -> \$N_OUT (\$N_COLLAPSED collapsed, \${PCT}%)"
    """
}
