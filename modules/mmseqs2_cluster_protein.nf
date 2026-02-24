/*
 * MMSEQS2_CLUSTER_PROTEIN — protein-level dedup at 95% amino acid identity
 * Catches near-identical proteins from different Corset clusters (recent
 * paralogs, allelic variants in polyploid Poaceae).  Directly reduces
 * duplicated BUSCOs.  Outputs a cluster map for optional tximport
 * aggregation of Salmon counts.
 */

process MMSEQS2_CLUSTER_PROTEIN {
    tag "${species_label}"

    input:
    path(protein_fasta)
    val(species_label)

    output:
    path("${species_label}.faa"),        emit: faa
    path("protein_cluster_map.tsv"),     emit: map
    path("protein_cluster_stats.txt"),   emit: stats

    script:
    """
    # Run MMseqs2 easy-cluster at 95% amino acid identity
    mmseqs easy-cluster \\
        ${protein_fasta} \\
        prot_clust \\
        tmp_prot_clust \\
        --min-seq-id 0.95 \\
        -c 0.8 \\
        --cov-mode 0 \\
        --threads ${task.cpus}

    # Representative proteins (deduped set)
    cp prot_clust_rep_seq.fasta ${species_label}.faa

    # Build cluster map: member_gene_id → representative_gene_id
    # (from prot_clust_cluster.tsv: rep_id \\t member_id)
    echo -e "gene_id\\tprotein_group_id" > protein_cluster_map.tsv
    awk -F'\\t' '{
        rep = \$1; sub(/ .*/, "", rep)
        mem = \$2; sub(/ .*/, "", mem)
        print mem "\\t" rep
    }' prot_clust_cluster.tsv >> protein_cluster_map.tsv

    # Report
    N_IN=\$(grep -c '^>' ${protein_fasta})
    N_OUT=\$(grep -c '^>' ${species_label}.faa)
    N_COLLAPSED=\$((N_IN - N_OUT))
    PCT=\$(awk "BEGIN{printf \\"%.1f\\", (\$N_IN > 0 ? \$N_COLLAPSED/\$N_IN*100 : 0)}")

    cat > protein_cluster_stats.txt <<EOF
MMseqs2 protein dedup (95% aa identity, 80% coverage)
Input proteins:      \$N_IN
Output proteins:     \$N_OUT
Collapsed:           \$N_COLLAPSED (\${PCT}%)
EOF

    echo "Protein dedup: \$N_IN -> \$N_OUT (\$N_COLLAPSED collapsed, \${PCT}%)"
    """
}
