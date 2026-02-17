/*
 * MMSEQS2_TAXONOMY — single-process taxonomy classification + filtering
 *
 * Copies TrEMBL DB to node-local SSD ($TMPDIR) to avoid NFS I/O bottleneck,
 * runs taxonomy on all SuperTranscripts at once, and cleans up local storage.
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
    # ── Use node-local SSD for all heavy I/O ──
    LOCAL=\${TMPDIR:-/tmp}/mmseqs2_taxonomy_\$\$
    mkdir -p \$LOCAL

    echo "Copying TrEMBL DB to local storage (\$LOCAL)..."
    for f in ${taxonomy_db}*; do
        cp "\$f" \$LOCAL/
    done
    LOCAL_DB=\$LOCAL/\$(basename ${taxonomy_db})
    echo "DB copy complete (\$(du -sh \$LOCAL | cut -f1))."

    # 1. Create MMseqs2 query DB (on local SSD)
    mmseqs createdb ${supertranscripts_fasta} \$LOCAL/queryDB

    # 2. Taxonomy assignment via LCA against UniProt/TrEMBL
    MEM_GB=\$(( ${task.memory.toGiga()} * 85 / 100 ))
    mmseqs taxonomy \$LOCAL/queryDB \$LOCAL_DB \$LOCAL/taxResult \$LOCAL/tmp_tax \\
        --tax-lineage 1 \\
        -s ${search_sens} \\
        --split-memory-limit \${MEM_GB}G \\
        --threads ${task.cpus}

    # 3. Export LCA results as TSV (write to work dir)
    mmseqs createtsv \$LOCAL/queryDB \$LOCAL/taxResult taxRes_lca.tsv

    # 4. Filter taxonomy: keep only target taxon and descendants
    mmseqs filtertaxdb \$LOCAL_DB \$LOCAL/taxResult \$LOCAL/filteredTaxResult \\
        --taxon-list ${filter_taxon}

    # 5. Extract matching query sequences
    mmseqs createsubdb \$LOCAL/filteredTaxResult \$LOCAL/queryDB \$LOCAL/filteredDB
    mmseqs createsubdb \$LOCAL/filteredTaxResult \$LOCAL/queryDB_h \$LOCAL/filteredDB_h

    # 6. Convert filtered DB back to FASTA (write to work dir)
    mmseqs convert2fasta \$LOCAL/filteredDB supertranscripts_filtered.fasta

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

    # 8. Cleanup local storage (also cleaned by SLURM epilog, but be tidy)
    rm -rf \$LOCAL
    """
}
