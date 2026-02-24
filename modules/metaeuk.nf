/*
 * METAEUK_PREDICT — homology-based gene prediction via MetaEuk easy-predict
 * Runs against SwissProt DB (MMseqs2 format).  Selects best protein per gene
 * (highest alignment score → lowest evalue → longest).
 */

process METAEUK_PREDICT {
    tag "${species_label}"

    input:
    path(fasta)
    val(swissprot_db)
    val(species_label)

    output:
    path("metaeuk_out.fas"),   emit: fas
    path("metaeuk_raw.gff"),   emit: gff

    script:
    """
    # Copy SwissProt DB to CWD (node-local SSD when scratch is enabled)
    for f in ${swissprot_db}*; do
        cp "\$f" .
    done
    LOCAL_DB=./\$(basename ${swissprot_db})

    # Run MetaEuk easy-predict
    metaeuk easy-predict \\
        ${fasta} \\
        \$LOCAL_DB \\
        metaeuk_out \\
        tmp_metaeuk \\
        --protein 1 \\
        --strand 2 \\
        -s 5.7 \\
        --min-length 30 \\
        --metaeuk-eval 0.001 \\
        --metaeuk-tcov 0.3 \\
        --max-intron 50000 \\
        --threads ${task.cpus} \\
        -v 2

    # Copy GFF before post-processing
    cp metaeuk_out.gff metaeuk_raw.gff
    """
}

/*
 * METAEUK_SELECT_BEST — pick best MetaEuk protein per gene
 * Runs in a Python-capable container (biopython) since the MetaEuk
 * biocontainer has no Python.
 */

process METAEUK_SELECT_BEST {
    tag "${species_label}"

    input:
    path(fas)
    val(species_label)

    output:
    path("metaeuk_best.faa"),  emit: faa
    path("metaeuk_map.tsv"),   emit: map

    script:
    """
    metaeuk_select_best.py ${fas} metaeuk_best.faa metaeuk_map.tsv
    """
}
