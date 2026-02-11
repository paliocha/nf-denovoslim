/*
 * Corset — hierarchical transcript-to-gene clustering with condition-aware
 *          paralog splitting, using Salmon equivalence classes
 */

process CORSET {
    label 'process_high'
    tag "${params.species_label}"

    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    path(quant_dirs)
    val(sample_conditions)   // list of maps: [[sample_id: X, condition: Y], ...]

    output:
    path("corset-clusters.txt"), emit: clust
    path("corset-counts.txt"),   emit: counts

    script:
    // Build -g (conditions) and -n (sample names) in the same order as sample_conditions
    def groups = sample_conditions.collect { it.condition }.join(',')
    def names  = sample_conditions.collect { it.sample_id }.join(',')
    def eq_files = sample_conditions.collect { "${it.sample_id}_quant/aux_info/eq_classes.txt" }.join(' ')
    """
    # Decompress eq_classes if Salmon produced gzipped output
    for d in *_quant/aux_info; do
        if [ -f "\$d/eq_classes.txt.gz" ] && [ ! -f "\$d/eq_classes.txt" ]; then
            gunzip -k "\$d/eq_classes.txt.gz"
        fi
    done

    corset \\
        -i salmon_eq_classes \\
        -g ${groups} \\
        -n ${names} \\
        -p corset \\
        ${eq_files}
    """
}
