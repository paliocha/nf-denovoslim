/*
 * PSAURON_GMST — score GeneMarkS-T CDS predictions with PSAURON
 * Uses nucleotide CDS input (no -p flag) for optimal accuracy.
 * PSAURON is bundled with TransDecoder2.
 */

process PSAURON_GMST {
    tag "${species_label}"

    input:
    path(gmst_fnn)
    val(species_label)

    output:
    path("gmst_psauron.csv"), emit: scores

    script:
    """
    # PSAURON scores CDS nucleotide sequences for coding potential
    # Input: nucleotide CDS FASTA (preferred over protein for accuracy)
    psauron -i ${gmst_fnn} --use-cpu -o gmst_psauron.csv

    N_SCORED=\$(tail -n +2 gmst_psauron.csv | grep -c ',' || echo 0)
    echo "PSAURON scored \$N_SCORED GeneMarkS-T CDS sequences"
    """
}
