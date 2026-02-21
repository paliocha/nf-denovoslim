/*
 * MMseqs2 taxonomy classification + filtering
 * DB is copied to CWD (node-local SSD when scratch is enabled).
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
    # Copy taxonomy DB to CWD (node-local SSD when scratch is enabled)
    for f in ${taxonomy_db}*; do
        cp "\$f" .
    done
    LOCAL_DB=./\$(basename ${taxonomy_db})

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

    # 6b. Rescue sequences with no database hit (taxon ID 0 = unclassified).
    #     Novel or species-specific plant genes lacking a UniRef90 match receive
    #     taxon ID 0 and are NOT kept by filtertaxdb, so they must be appended
    #     explicitly.  Contamination sequences (e.g. bacteria, fungi) will have
    #     a real taxonomy assignment and are correctly excluded by step 4.
    awk -F'\t' '\$2 == 0 { print \$1 }' taxRes_lca.tsv > no_hit_ids.txt
    if [ -s no_hit_ids.txt ]; then
        awk 'NR==FNR { ids[\$1]=1; next }
             /^>/ { p = (substr(\$1,2) in ids) }
             p' no_hit_ids.txt ${supertranscripts_fasta} >> supertranscripts_filtered.fasta
    fi

    # 7. Statistics
    TOTAL=\$(grep -c '^>' ${supertranscripts_fasta})
    KEPT=\$(grep -c '^>' supertranscripts_filtered.fasta || echo 0)
    NOHIT=\$(wc -l < no_hit_ids.txt)
    KEPT_TAXON=\$((KEPT - NOHIT))
    REMOVED=\$((TOTAL - KEPT))

    cat > taxonomy_filter_stats.txt <<EOF
Taxonomy filter: NCBI taxon ID ${filter_taxon}
Total SuperTranscripts:         \$TOTAL
Kept (total):                   \$KEPT (\$(awk "BEGIN{printf \\"%.1f\\", \$KEPT/\$TOTAL*100}")%)
  - Streptophyta-assigned:      \$KEPT_TAXON
  - No taxonomy hit (rescued):  \$NOHIT
Removed (non-plant):            \$REMOVED (\$(awk "BEGIN{printf \\"%.1f\\", \$REMOVED/\$TOTAL*100}")%)
EOF

    echo "Taxonomy filter: kept \$KEPT / \$TOTAL SuperTranscripts (\$NOHIT with no hit rescued)"
    """
}
