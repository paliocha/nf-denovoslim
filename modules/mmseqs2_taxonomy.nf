/*
 * MMseqs2 taxonomy classification + filtering
 * DB is copied to CWD (node-local SSD when scratch is enabled).
 *
 * Filter logic:
 *   KEEP  — descendants of ${filter_taxon} (Viridiplantae by default)
 *   KEEP  — taxon ID 0 (no database hit: novel / species-specific genes)
 *   REMOVE — everything else (Bacteria, Archaea, Fungi, Metazoa, Viruses,
 *            non-plant Eukaryota such as Oomycetes/Stramenopiles)
 *
 * CRITICAL: mmseqs filtertaxdb does NOT delete non-matching entries; it
 * zeros their data (size → 1 byte) but keeps the key.  A plain
 * `createsubdb filteredTaxResult queryDB` would therefore extract ALL
 * sequences.  We use `awk '$3 != 1'` on the .index to collect only
 * the keys whose data was NOT zeroed, then pass that ID list to
 * createsubdb.  (See MMseqs2 wiki "filtertaxdb" example.)
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
    path("taxonomy_breakdown.tsv"),          emit: breakdown

    script:
    """
    # Copy taxonomy DB to CWD (node-local SSD when scratch is enabled)
    for f in ${taxonomy_db}*; do
        cp "\$f" .
    done
    LOCAL_DB=./\$(basename ${taxonomy_db})

    # 1. Create MMseqs2 query DB
    mmseqs createdb ${supertranscripts_fasta} queryDB

    # 2. Taxonomy assignment via LCA
    MEM_GB=\$(( ${task.memory.toGiga()} * 85 / 100 ))
    mmseqs taxonomy queryDB \$LOCAL_DB taxResult tmp_tax \\
        --tax-lineage 1 \\
        -s ${search_sens} \\
        --split-memory-limit \${MEM_GB}G \\
        --threads ${task.cpus}

    # 3. Export LCA results as TSV (with lineage in column 9)
    mmseqs createtsv queryDB taxResult taxRes_lca.tsv

    # 4. Taxonomy filter: keep Viridiplantae + no-hit only
    #    ${filter_taxon} = 33090 (Viridiplantae) by default
    #    0               = sequences with no taxonomy assignment
    mmseqs filtertaxdb \$LOCAL_DB taxResult filteredTaxResult \\
        --taxon-list ${filter_taxon},0

    # 5. Extract ONLY keys whose data was NOT zeroed by filtertaxdb
    #    (column 3 of the .index = data size; size 1 = zeroed/removed)
    awk '\$3 != 1 { print \$1 }' filteredTaxResult.index > kept_ids.list

    # 6. Build filtered sequence DB from the kept-ID list
    mmseqs createsubdb kept_ids.list queryDB   filteredDB
    mmseqs createsubdb kept_ids.list queryDB_h filteredDB_h

    # 7. Convert filtered DB back to FASTA
    mmseqs convert2fasta filteredDB supertranscripts_filtered.fasta

    # 8. Taxonomy lineage breakdown table
    #    Categorises every SuperTranscript by broad taxonomic group using
    #    the lineage string (column 9) from taxRes_lca.tsv, then marks
    #    each category as KEPT or REMOVED based on the filter.
    awk -F'\\t' '
    {
        taxid = \$2
        lineage = \$9

        if (taxid == 0)
            cat = "No taxonomy hit"
        else if (lineage ~ /Viridiplantae/) {
            if (lineage ~ /Streptophyta/)
                cat = "Streptophyta"
            else
                cat = "Viridiplantae (non-Streptophyta)"
        }
        else if (lineage ~ /Fungi/)
            cat = "Fungi"
        else if (lineage ~ /Metazoa/)
            cat = "Metazoa"
        else if (lineage ~ /Bacteria/)
            cat = "Bacteria"
        else if (lineage ~ /Archaea/)
            cat = "Archaea"
        else if (lineage ~ /Virus/)
            cat = "Viruses"
        else if (lineage ~ /Eukaryota/)
            cat = "Other Eukaryota"
        else
            cat = "Other"

        counts[cat]++
    }
    END {
        OFS = "\\t"
        print "Category", "Count", "Percent", "Status"

        total = 0
        for (c in counts) total += counts[c]

        n = split("Streptophyta,Viridiplantae (non-Streptophyta),No taxonomy hit,Other Eukaryota,Fungi,Metazoa,Bacteria,Archaea,Viruses,Other", order, ",")
        for (i = 1; i <= n; i++) {
            c = order[i]
            if (c in counts) {
                pct = sprintf("%.1f", counts[c] / total * 100)
                if (c == "Streptophyta" || c == "Viridiplantae (non-Streptophyta)" || c == "No taxonomy hit")
                    status = "KEPT"
                else
                    status = "REMOVED"
                print c, counts[c], pct "%", status
            }
        }
        print "TOTAL", total, "100.0%", "-"
    }
    ' taxRes_lca.tsv > taxonomy_breakdown.tsv

    # 9. Summary statistics
    TOTAL=\$(grep -c '^>' ${supertranscripts_fasta})
    KEPT=\$(grep -c '^>' supertranscripts_filtered.fasta || echo 0)
    REMOVED=\$((TOTAL - KEPT))

    # Counts from breakdown
    NOHIT=\$(awk -F'\\t' '\$1 == "No taxonomy hit" { print \$2 }' taxonomy_breakdown.tsv)
    NOHIT=\${NOHIT:-0}
    PLANT=\$((KEPT - NOHIT))

    cat > taxonomy_filter_stats.txt <<EOF
Taxonomy filter: ${filter_taxon} (Viridiplantae) + 0 (no-hit)
Total SuperTranscripts:         \$TOTAL
Kept (total):                   \$KEPT (\$(awk "BEGIN{printf \\"%.1f\\", \$KEPT/\$TOTAL*100}")%)
  - Viridiplantae-assigned:     \$PLANT
  - No taxonomy hit (rescued):  \$NOHIT
Removed (non-Viridiplantae):    \$REMOVED (\$(awk "BEGIN{printf \\"%.1f\\", \$REMOVED/\$TOTAL*100}")%)

--- Lineage Breakdown ---
\$(cat taxonomy_breakdown.tsv | column -t -s '\t')
EOF

    echo "Taxonomy filter: kept \$KEPT / \$TOTAL SuperTranscripts (\$REMOVED removed, \$NOHIT no-hit rescued)"
    """
}
