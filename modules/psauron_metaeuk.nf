/*
 * PSAURON_METAEUK — re-score MetaEuk protein predictions with PSAURON
 * Enables quality-based filtering in the merge step and fair comparison
 * with TD2 predictions.  PSAURON is bundled with TransDecoder2.
 */

process PSAURON_METAEUK {
    tag "${species_label}"

    input:
    path(metaeuk_faa)
    val(species_label)

    output:
    path("metaeuk_psauron.csv"), emit: scores

    script:
    """
    # PSAURON scores protein sequences for coding potential
    # Input: protein FASTA (-p flag); Output: CSV with description,in-frame_score columns
    psauron -p -i ${metaeuk_faa} --use-cpu -o metaeuk_psauron.csv

    N_SCORED=\$(tail -n +2 metaeuk_psauron.csv | grep -c ',' || echo 0)
    echo "PSAURON scored \$N_SCORED MetaEuk proteins"
    """
}
