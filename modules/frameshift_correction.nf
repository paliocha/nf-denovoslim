/*
 * FRAMESHIFT_CORRECTION — detect and correct assembly frameshifts
 *
 * Uses Diamond blastx with -F 15 (frameshift-tolerant alignment) to identify
 * frameshifts introduced during Trinity assembly. The BTOP string from Diamond
 * encodes frameshift positions, which are then corrected in the nucleotide
 * sequence by a custom Python script.
 *
 * Based on: Leder et al. (2021) J Evol Biol 34:138 — Diamond BTOP parsing
 * for frameshift removal in de novo transcriptomes.
 *
 * Input:  taxonomy-filtered SuperTranscripts FASTA
 * Output: corrected SuperTranscripts FASTA (frameshifts fixed in-place)
 *
 * Requires params.diamond_db — path to a pre-built Diamond .dmnd database
 * or a protein FASTA file (will be indexed automatically).
 *
 * Build a Diamond DB from SwissProt:
 *   wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
 *   diamond makedb --in uniprot_sprot.fasta.gz -d swissprot
 */

process FRAMESHIFT_CORRECTION {
    label 'process_medium'
    tag "${params.species_label}"

    publishDir "${params.outdir}/frameshift_correction", mode: 'copy',
        pattern: 'frameshift_stats.txt'

    input:
    path(supertranscripts_fasta)

    output:
    path("supertranscripts_corrected.fasta"), emit: fasta
    path("frameshift_stats.txt"),             emit: stats

    script:
    """
    # 1. Prepare Diamond DB
    #    If .dmnd → use directly; if FASTA → build DB; if null → error
    DB_PATH="${params.diamond_db}"
    if [ -z "\$DB_PATH" ] || [ "\$DB_PATH" = "null" ]; then
        echo "ERROR: --diamond_db is required. Provide a .dmnd file or protein FASTA." >&2
        echo "Build from SwissProt: diamond makedb --in uniprot_sprot.fasta.gz -d swissprot" >&2
        exit 1
    fi

    if [[ "\$DB_PATH" == *.dmnd ]]; then
        DIAMOND_DB="\$DB_PATH"
    else
        # Assume protein FASTA — build Diamond DB
        diamond makedb --in \$DB_PATH -d diamond_ref -p ${task.cpus}
        DIAMOND_DB="diamond_ref"
    fi

    # 2. Run Diamond blastx with frameshift-tolerant alignment
    diamond blastx \\
        -F 15 \\
        --sensitive \\
        --top 1 \\
        --min-score 50 \\
        -d \$DIAMOND_DB \\
        -q ${supertranscripts_fasta} \\
        --outfmt 6 qseqid qstart qend qlen qframe btop \\
        -p ${task.cpus} \\
        -o diamond_fs.tsv

    # 3. Parse BTOP strings and correct frameshifts
    correct_frameshifts.py \\
        ${supertranscripts_fasta} \\
        diamond_fs.tsv \\
        supertranscripts_corrected.fasta \\
        | tee frameshift_stats.txt
    """
}
