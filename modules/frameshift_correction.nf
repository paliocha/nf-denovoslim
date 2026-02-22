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

    # Run Diamond blastx with frameshift-tolerant alignment
    # --iterate is incompatible with -F (frameshift): iterate's non-linear steps
    # (default, sensitive) use full matrix extension which -F does not support.
    diamond blastx \\
        -F 15 \\
        --sensitive \\
        --strand plus \\
        --top 1 \\
        --min-score 50 \\
        -b 4 -c 1 \\
        -d \$LOCAL_DB \\
        -q ${supertranscripts_fasta} \\
        --outfmt 6 qseqid qstart qend qlen qframe btop \\
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
