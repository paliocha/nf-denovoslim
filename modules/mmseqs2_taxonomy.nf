/*
 * MMSEQS2_TAXONOMY — single-process taxonomy classification + filtering
 *
 * With scratch enabled (Orion), the task CWD is already on node-local SSD.
 * The taxonomy DB is copied to CWD for fast memory-mapped access.
 * With 1.5 TB RAM per node and 3.5 TB local SSD, this is faster than chunking
 * because the DB is loaded only once instead of N times.
 */

process MMSEQS2_TAXONOMY {
    tag "${species_label}"

    input:
    path(supertranscripts_fasta)
    val(taxonomy_db)
    val(search_sens)
    val(filter_taxon)
    val(species_label)

    output:
    path("supertranscripts_filtered.fasta"), emit: fasta
    path("taxRes_lca.tsv"),                  emit: lca_tsv
    path("taxonomy_filter_stats.txt"),       emit: stats

    script:
    """
    # ── Copy taxonomy DB to CWD (node-local SSD when scratch is enabled) ──
    echo "Copying taxonomy DB to local storage..."
    for f in ${taxonomy_db}*; do
        cp "\$f" .
    done
    LOCAL_DB=./\$(basename ${taxonomy_db})
    echo "DB copy complete (\$(du -sh \$LOCAL_DB* | tail -1 | cut -f1))."

    # 1. Create MMseqs2 query DB
    mmseqs createdb ${supertranscripts_fasta} queryDB

    # 2. Taxonomy assignment via LCA against UniRef90
    MEM_GB=\$(( ${task.memory.toGiga()} * 85 / 100 ))
    mmseqs taxonomy queryDB \$LOCAL_DB taxResult tmp_tax \\
        --tax-lineage 1 \\
        -s ${search_sens} \\
        --split-memory-limit \${MEM_GB}G \\
        --threads ${task.cpus}

    # 3. Export LCA results as TSV
    mmseqs createtsv queryDB taxResult taxRes_lca.tsv

    # 4. Filter taxonomy: keep only target taxon and descendants
    mmseqs filtertaxdb \$LOCAL_DB taxResult filteredTaxResult \\
        --taxon-list ${filter_taxon}

    # 5. Extract matching query sequences
    mmseqs createsubdb filteredTaxResult queryDB filteredDB
    mmseqs createsubdb filteredTaxResult queryDB_h filteredDB_h

    # 6. Convert filtered DB back to FASTA
    mmseqs convert2fasta filteredDB supertranscripts_filtered.fasta

    # 7. Statistics
    TOTAL=\$(grep -c '^>' ${supertranscripts_fasta})
    KEPT=\$(grep -c '^>' supertranscripts_filtered.fasta || echo 0)
    REMOVED=\$((TOTAL - KEPT))

    cat > taxonomy_filter_stats.txt <<EOF
Taxonomy filter: NCBI taxon ID ${filter_taxon}
Total SuperTranscripts:    \$TOTAL
Matching (kept):           \$KEPT (\$(awk "BEGIN{printf \\"%.1f\\", \$KEPT/\$TOTAL*100}")%)
Non-matching (removed):    \$REMOVED (\$(awk "BEGIN{printf \\"%.1f\\", \$REMOVED/\$TOTAL*100}")%)
EOF

    echo "Taxonomy filter: kept \$KEPT / \$TOTAL SuperTranscripts (removed \$REMOVED)"
    """
}
