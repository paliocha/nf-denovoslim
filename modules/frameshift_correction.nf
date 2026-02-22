/*
 * Frameshift correction — Diamond blastx + Python correction
 * Diamond DB is copied to CWD (node-local SSD when scratch is enabled).
 * Python correction is a separate process so Diamond results stay cached.
 */

process DIAMOND_BLASTX {
    tag "${species_label}"

    input:
    path(supertranscripts_fasta)
    val(diamond_db)
    val(species_label)

    output:
    path("diamond_frameshift.tsv"), emit: tsv

    script:
    """
    # Copy Diamond DB to CWD (node-local SSD when scratch is enabled)
    cp "${diamond_db}" .
    LOCAL_DB=./\$(basename "${diamond_db}")

    # --- Dynamically compute optimal -b and -c from available RAM ----------
    # Diamond memory model: ~6b + 6b/c GB  (b = block size in billions of
    # letters, c = index chunks).  Goal: fit the entire DB in one block
    # (c=1) for minimum passes.  If RAM is tight, increase c to reduce
    # per-pass index memory at the cost of more passes.
    MEM_GB=${task.memory.toGiga()}
    DB_LETTERS=\$(diamond dbinfo -d \$LOCAL_DB 2>/dev/null | awk '/^letters/{print \$2}')
    # Minimum -b to hold the whole DB in one reference block
    MIN_B=\$(( (DB_LETTERS + 999999999) / 1000000000 ))   # ceil(letters / 1e9)

    # Available for Diamond: 85% of allocated RAM (leave headroom for OS + staging)
    AVAIL_GB=\$(( MEM_GB * 85 / 100 ))

    # With c=1: need ~12*b GB.  With c=N: need ~6*b*(1 + 1/N) ≈ 6*b for large N.
    if (( AVAIL_GB >= 12 * MIN_B )); then
        BLOCK=\$MIN_B
        CHUNKS=1           # 1 block × 1 chunk × 16 shapes = 16 passes (minimum)
    elif (( AVAIL_GB >= 8 * MIN_B )); then
        BLOCK=\$MIN_B
        CHUNKS=\$(( (6 * MIN_B + AVAIL_GB - 6 * MIN_B - 1) / (AVAIL_GB - 6 * MIN_B) ))
        (( CHUNKS < 2 )) && CHUNKS=2
    else
        # Not enough RAM for 1 block — let Diamond auto-split
        BLOCK=\$(( AVAIL_GB / 12 ))
        (( BLOCK < 1 )) && BLOCK=1
        CHUNKS=1
    fi

    echo "Diamond tuning: DB_LETTERS=\${DB_LETTERS}, MIN_B=\${MIN_B}, MEM_GB=\${MEM_GB}, AVAIL_GB=\${AVAIL_GB} → -b \${BLOCK} -c \${CHUNKS}"

    # Run Diamond blastx with frameshift-tolerant alignment.
    # --iterate is incompatible with -F: non-linear steps use full matrix
    # extension which frameshift mode does not support.
    diamond blastx \\
        -F 15 \\
        --sensitive \\
        --strand plus \\
        --top 5 \\
        --min-score 60 \\
        --id 60 \\
        --unal 1 \\
        -b \$BLOCK -c \$CHUNKS \\
        -d \$LOCAL_DB \\
        -q ${supertranscripts_fasta} \\
        --outfmt 6 qseqid qstart qend qlen qframe btop \\
            sseqid slen sstart send evalue bitscore score length pident qseq \\
        -p ${task.cpus} \\
        --tmpdir . \\
        -o diamond_frameshift.tsv
    """
}

process CORRECT_FRAMESHIFTS {
    tag "${species_label}"

    input:
    path(supertranscripts_fasta)
    path(diamond_tsv)
    val(species_label)

    output:
    path("supertranscripts_corrected.fasta"), emit: fasta
    path("frameshift_stats.txt"),             emit: stats

    script:
    """
    correct_frameshifts.py \\
        ${supertranscripts_fasta} \\
        ${diamond_tsv} \\
        supertranscripts_corrected.fasta \\
        > frameshift_stats.txt
    """
}
