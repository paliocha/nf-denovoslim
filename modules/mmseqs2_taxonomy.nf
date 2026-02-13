/*
 * MMSEQS2_TAXONOMY — classify and filter SuperTranscripts by taxonomy
 *
 * Uses MMseqs2's native taxonomy workflow:
 *   1. createdb    — convert FASTA to MMseqs2 DB
 *   2. taxonomy    — LCA assignment via 6-frame search against UniProt/TrEMBL
 *   3. filtertaxdb — keep only hits under target taxon (NCBI taxonomy hierarchy)
 *   4. createsubdb — extract matching sequences
 *   5. convert2fasta — back to FASTA
 *
 * This is more robust than grepping lineage strings because filtertaxdb
 * uses the NCBI taxonomy tree to include all descendants of the target taxon.
 */

process MMSEQS2_TAXONOMY {
    label 'process_high'
    tag "${params.species_label}"

    publishDir "${params.outdir}/taxonomy", mode: 'copy'

    input:
    path(supertranscripts_fasta)

    output:
    path("supertranscripts_filtered.fasta"), emit: fasta
    path("taxRes_lca.tsv"),                  emit: lca_tsv
    path("taxRes_report"),                   emit: report
    path("taxonomy_filter_stats.txt"),       emit: stats

    script:
    """
    # 1. Create MMseqs2 query DB from SuperTranscripts
    mmseqs createdb ${supertranscripts_fasta} queryDB

    # 2. Taxonomy assignment via LCA against UniProt/TrEMBL (6-frame translation)
    MEM_GB=\$(echo "${task.memory.toGiga()} * 85 / 100" | bc)
    mmseqs taxonomy queryDB ${params.mmseqs2_taxonomy_db} taxResult tmp_tax \\
        --tax-lineage 1 \\
        -s ${params.mmseqs2_search_sens} \\
        --split-memory-limit \${MEM_GB}G \\
        --threads ${task.cpus}

    # 3. Export full LCA report as TSV (for downstream inspection)
    mmseqs createtsv queryDB taxResult taxRes_lca.tsv

    # 4. Generate Krona-compatible taxonomy report
    mmseqs taxonomyreport ${params.mmseqs2_taxonomy_db} taxResult taxRes_report

    # 5. Filter taxonomy results: keep only target taxon and all descendants
    mmseqs filtertaxdb ${params.mmseqs2_taxonomy_db} taxResult filteredTaxResult \\
        --taxon-list ${params.filter_taxon}

    # 6. Extract matching query sequences
    mmseqs createsubdb filteredTaxResult queryDB filteredDB
    mmseqs createsubdb filteredTaxResult queryDB_h filteredDB_h

    # 7. Convert filtered DB back to FASTA
    mmseqs convert2fasta filteredDB supertranscripts_filtered.fasta

    # 8. Compute filter statistics
    TOTAL=\$(grep -c '^>' ${supertranscripts_fasta})
    KEPT=\$(grep -c '^>' supertranscripts_filtered.fasta || echo 0)
    REMOVED=\$((TOTAL - KEPT))

    cat > taxonomy_filter_stats.txt <<EOF
Taxonomy filter: NCBI taxon ID ${params.filter_taxon}
Total SuperTranscripts:    \$TOTAL
Matching (kept):           \$KEPT (\$(awk "BEGIN{printf \\"%.1f\\", \$KEPT/\$TOTAL*100}")%)
Non-matching (removed):    \$REMOVED (\$(awk "BEGIN{printf \\"%.1f\\", \$REMOVED/\$TOTAL*100}")%)
EOF

    echo "Taxonomy filter: kept \$KEPT / \$TOTAL SuperTranscripts (removed \$REMOVED non-taxid-${params.filter_taxon})"
    """
}
