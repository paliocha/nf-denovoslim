/*
 * FRAMESHIFT_CORRECTION — single-process Diamond blastx + Python correction
 *
 * Diamond blastx copies uniref90.dmnd to node-local SSD ($TMPDIR) to avoid
 * NFS memory-mapping issues (SIGBUS). Python correction is a separate process
 * so Diamond results stay cached if the lightweight correction step needs re-running.
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
    # ── Use node-local SSD to avoid NFS SIGBUS on memory-mapped .dmnd ──
    LOCAL=\${TMPDIR:-/tmp}/diamond_blastx_\$\$
    mkdir -p \$LOCAL

    # Copy DB to local SSD (~86 GB, ~1-2 min)
    echo "Copying Diamond DB to local storage (\$LOCAL)..."
    cp "${diamond_db}" \$LOCAL/
    LOCAL_DB=\$LOCAL/\$(basename "${diamond_db}")
    echo "DB copy complete (\$(du -sh \$LOCAL | cut -f1))."

    # Run Diamond blastx with frameshift-tolerant alignment
    diamond blastx \\
        -F 15 \\
        --sensitive \\
        --top 1 \\
        --min-score 50 \\
        -d \$LOCAL_DB \\
        -q ${supertranscripts_fasta} \\
        --outfmt 6 qseqid qstart qend qlen qframe btop \\
        -p ${task.cpus} \\
        --tmpdir \$LOCAL \\
        -o diamond_frameshift.tsv

    echo "Diamond blastx complete: \$(wc -l < diamond_frameshift.tsv) alignments"

    # Cleanup local storage (also cleaned by SLURM epilog)
    rm -rf \$LOCAL
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

    echo "Frameshift correction: \$(grep -c '^>' supertranscripts_corrected.fasta) sequences"
    """
}
