/*
 * MMSEQS2_TAXONOMY with chunking for parallel execution
 *
 * Splits SuperTranscripts into chunks, runs taxonomy classification and filtering
 * in parallel, then merges results. Provides 5-10× speedup (24h → 3-4h).
 */

process SPLIT_SUPERTRANSCRIPTS {
    tag "${params.species_label}"

    input:
    path(supertranscripts_fasta)

    output:
    path("chunks/st_chunk_*.fasta"), emit: chunks

    script:
    def chunk_size = params.taxonomy_chunk_size ?: 2000
    """
    mkdir -p chunks

    # Split into chunks of ${chunk_size} sequences each
    awk -v size=${chunk_size} -v pre=chunks/st_chunk_ '
        /^>/ {
            if (n % size == 0) {
                file = sprintf("%s%03d.fasta", pre, int(n / size))
            }
            n++
        }
        {print > file}
    ' ${supertranscripts_fasta}

    # Ensure at least one chunk exists even for small inputs
    if [ ! -f chunks/st_chunk_000.fasta ]; then
        cp ${supertranscripts_fasta} chunks/st_chunk_000.fasta
    fi

    echo "Split \$(grep -c '^>' ${supertranscripts_fasta}) SuperTranscripts into \$(ls chunks/*.fasta | wc -l) chunks"
    """
}

process MMSEQS2_TAXONOMY_CHUNK {
    tag "${params.species_label}_chunk_${chunk_idx}"

    input:
    tuple val(chunk_idx), path(chunk_fasta)

    output:
    tuple val(chunk_idx), path("filtered_${chunk_idx}.fasta"), emit: fasta
    tuple val(chunk_idx), path("tax_${chunk_idx}.tsv"),        emit: tsv
    path("stats_${chunk_idx}.txt"),                            emit: stats

    script:
    """
    # 1. Create MMseqs2 query DB from chunk
    mmseqs createdb ${chunk_fasta} queryDB

    # 2. Taxonomy assignment via LCA against UniProt/TrEMBL
    MEM_GB=\$(( ${task.memory.toGiga()} * 85 / 100 ))
    mmseqs taxonomy queryDB ${params.mmseqs2_taxonomy_db} taxResult tmp_tax \\
        --tax-lineage 1 \\
        -s ${params.mmseqs2_search_sens} \\
        --split-memory-limit \${MEM_GB}G \\
        --threads ${task.cpus}

    # 3. Export LCA results as TSV
    mmseqs createtsv queryDB taxResult tax_${chunk_idx}.tsv

    # 4. Filter taxonomy results: keep only target taxon and descendants
    mmseqs filtertaxdb ${params.mmseqs2_taxonomy_db} taxResult filteredTaxResult \\
        --taxon-list ${params.filter_taxon}

    # 5. Extract matching query sequences
    mmseqs createsubdb filteredTaxResult queryDB filteredDB
    mmseqs createsubdb filteredTaxResult queryDB_h filteredDB_h

    # 6. Convert filtered DB back to FASTA
    mmseqs convert2fasta filteredDB filtered_${chunk_idx}.fasta

    # 7. Chunk statistics
    TOTAL=\$(grep -c '^>' ${chunk_fasta})
    KEPT=\$(grep -c '^>' filtered_${chunk_idx}.fasta || echo 0)
    REMOVED=\$((TOTAL - KEPT))

    cat > stats_${chunk_idx}.txt <<EOF
Chunk ${chunk_idx}:
Total:   \$TOTAL
Kept:    \$KEPT
Removed: \$REMOVED
EOF

    echo "Chunk ${chunk_idx}: kept \$KEPT / \$TOTAL SuperTranscripts"
    """
}

process MERGE_TAXONOMY_RESULTS {
    tag "${params.species_label}"
    publishDir "${params.outdir}/taxonomy", mode: 'copy'

    input:
    path("filtered_*.fasta")
    path("tax_*.tsv")
    path("stats_*.txt")
    path(original_fasta)

    output:
    path("supertranscripts_filtered.fasta"), emit: fasta
    path("taxRes_lca.tsv"),                  emit: lca_tsv
    path("taxonomy_filter_stats.txt"),       emit: stats

    script:
    """
    # Merge filtered FASTA results
    cat filtered_*.fasta > supertranscripts_filtered.fasta

    # Merge taxonomy TSV results (sort by query ID for consistency)
    cat tax_*.tsv | sort -k1,1 > taxRes_lca.tsv

    # Compute combined filter statistics
    TOTAL=\$(grep -c '^>' ${original_fasta})
    KEPT=\$(grep -c '^>' supertranscripts_filtered.fasta || echo 0)
    REMOVED=\$((TOTAL - KEPT))

    cat > taxonomy_filter_stats.txt <<EOF
Taxonomy filter: NCBI taxon ID ${params.filter_taxon} (chunked processing)
Total SuperTranscripts:    \$TOTAL
Matching (kept):           \$KEPT (\$(awk "BEGIN{printf \\"%.1f\\", \$KEPT/\$TOTAL*100}")%)
Non-matching (removed):    \$REMOVED (\$(awk "BEGIN{printf \\"%.1f\\", \$REMOVED/\$TOTAL*100}")%)

Per-chunk statistics:
EOF

    # Append per-chunk stats
    cat stats_*.txt >> taxonomy_filter_stats.txt

    echo "Taxonomy filter (chunked): kept \$KEPT / \$TOTAL SuperTranscripts (removed \$REMOVED)"
    echo "Merged \$(ls filtered_*.fasta | wc -l) chunks"
    """
}
