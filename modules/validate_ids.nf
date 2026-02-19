/*
 * VALIDATE_IDS — check .faa IDs are a subset of Salmon gene IDs
 */

process VALIDATE_IDS {
    tag "${species_label}"

    input:
    path(quant_sf)
    path(faa)
    val(species_label)

    output:
    path("id_validation.txt"), emit: report

    script:
    """
    # Extract Salmon gene IDs
    cut -f1 ${quant_sf} | tail -n+2 | sort > ids_salmon.txt

    # Extract .faa protein IDs (first word only, strip description)
    grep "^>" ${faa} | awk '{sub(/^>/, ""); print \$1}' | sort > ids_faa.txt

    # Proteins must be a subset of Salmon genes
    ORPHANS=\$(comm -23 ids_faa.txt ids_salmon.txt | wc -l)
    if [ "\$ORPHANS" -gt 0 ]; then
        echo "ERROR: \${ORPHANS} protein IDs not found in Salmon quant!" >&2
        comm -23 ids_faa.txt ids_salmon.txt >&2
        exit 1
    fi

    # Report
    NO_ORF=\$(comm -23 ids_salmon.txt ids_faa.txt | wc -l)
    {
        echo "Genes with protein: \$(wc -l < ids_faa.txt)"
        echo "Genes without ORF (non-coding/fragmented): \${NO_ORF}"
        echo "Total SuperTranscript genes: \$(wc -l < ids_salmon.txt)"
    } | tee id_validation.txt
    """
}
