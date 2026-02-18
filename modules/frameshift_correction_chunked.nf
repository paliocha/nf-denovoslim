/*
 * FRAMESHIFT_CORRECTION with chunking for parallel execution
 *
 * Splits SuperTranscripts into chunks, runs Diamond blastx and frameshift
 * correction in parallel, then merges results. Provides 3-5x speedup
 * (8h -> ~2h) with better fault tolerance.
 *
 * Diamond and Python correction are separate processes so that Diamond
 * results stay cached if the Python step fails.
 */

process SPLIT_SUPERTRANSCRIPTS_FS {
    tag "${species_label}"

    input:
    path(supertranscripts_fasta)
    val(species_label)

    output:
    path("chunks/fs_chunk_*.fasta"), emit: chunks

    script:
    // Determine chunk count dynamically based on input size
    """
    mkdir -p chunks

    # Count total sequences
    TOTAL=\$(grep -c '^>' ${supertranscripts_fasta})

    # Aim for ~5 chunks (adjustable)
    CHUNK_SIZE=\$(( (\$TOTAL / 5) + 1 ))

    # Split into chunks
    awk -v size=\$CHUNK_SIZE -v pre=chunks/fs_chunk_ '
        /^>/ {
            if (n % size == 0) {
                file = sprintf("%s%02d.fasta", pre, int(n / size))
            }
            n++
        }
        {print > file}
    ' ${supertranscripts_fasta}

    # Ensure at least one chunk exists
    if [ ! -f chunks/fs_chunk_00.fasta ]; then
        cp ${supertranscripts_fasta} chunks/fs_chunk_00.fasta
    fi

    echo "Split \$TOTAL SuperTranscripts into \$(ls chunks/*.fasta | wc -l) chunks"
    """
}

process DIAMOND_BLASTX_CHUNK {
    tag "${species_label}_chunk_${chunk_idx}"

    input:
    tuple val(chunk_idx), path(chunk_fasta)
    val(diamond_db)
    val(species_label)

    output:
    tuple val(chunk_idx), path(chunk_fasta), path("diamond_fs_${chunk_idx}.tsv"), emit: results

    script:
    """
    # Prepare Diamond DB
    DB_PATH="${diamond_db}"
    if [ -z "\$DB_PATH" ] || [ "\$DB_PATH" = "null" ]; then
        echo "ERROR: --diamond_db is required." >&2
        exit 1
    fi

    if [[ "\$DB_PATH" == *.dmnd ]]; then
        DIAMOND_DB="\$DB_PATH"
    else
        diamond makedb --in \$DB_PATH -d diamond_ref -p ${task.cpus}
        DIAMOND_DB="diamond_ref"
    fi

    # Run Diamond blastx with frameshift-tolerant alignment
    diamond blastx \\
        -F 15 \\
        --sensitive \\
        --top 1 \\
        --min-score 50 \\
        -d \$DIAMOND_DB \\
        -q ${chunk_fasta} \\
        --outfmt 6 qseqid qstart qend qlen qframe btop \\
        -p ${task.cpus} \\
        -o diamond_fs_${chunk_idx}.tsv
    """
}

process CORRECT_FRAMESHIFTS_CHUNK {
    tag "${species_label}_chunk_${chunk_idx}"

    input:
    tuple val(chunk_idx), path(chunk_fasta), path(diamond_tsv)
    val(species_label)

    output:
    tuple val(chunk_idx), path("corrected_${chunk_idx}.fasta"), emit: fasta
    path("stats_${chunk_idx}.txt"),                              emit: stats

    script:
    """
    correct_frameshifts.py \\
        ${chunk_fasta} \\
        ${diamond_tsv} \\
        corrected_${chunk_idx}.fasta \\
        > stats_${chunk_idx}.txt

    echo "Chunk ${chunk_idx}: \$(grep -c '^>' corrected_${chunk_idx}.fasta || echo 0) sequences processed"
    """
}

process MERGE_FRAMESHIFT_RESULTS {
    tag "${species_label}"

    input:
    path("corrected_*.fasta")
    path("stats_*.txt")
    val(species_label)

    output:
    path("supertranscripts_corrected.fasta"), emit: fasta
    path("frameshift_stats.txt"),             emit: stats

    script:
    """
    # Merge corrected FASTA files
    cat corrected_*.fasta > supertranscripts_corrected.fasta

    # Merge statistics
    cat > frameshift_stats.txt <<EOF
Frameshift Correction Summary (chunked processing)
Total sequences: \$(grep -c '^>' supertranscripts_corrected.fasta)

Per-chunk statistics:
EOF

    # Append per-chunk stats
    cat stats_*.txt >> frameshift_stats.txt

    echo "Merged \$(ls corrected_*.fasta | wc -l) chunks into supertranscripts_corrected.fasta"
    """
}
