/*
 * Grouper — expression-aware transcript-to-gene clustering
 */

process GROUPER {
    label 'process_high'
    tag "${params.species_label}"

    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    path(quant_dirs)
    val(sample_conditions)   // list of maps: [[sample_id: X, condition: Y], ...]

    output:
    path("grouper_out/mag.flat.clust"), emit: clust
    path("grouper_out"),                emit: outdir

    script:
    // Build YAML config dynamically from sample metadata
    def conditions_set = sample_conditions.collect { it.condition }.unique().sort()
    def yaml_conditions = conditions_set.collect { "    - ${it}" }.join('\n')
    def yaml_samples = conditions_set.collect { cond ->
        def samples = sample_conditions
            .findAll { it.condition == cond }
            .collect { "        - ${it.sample_id}_quant" }
            .join('\n')
        "    ${cond}:\n${samples}"
    }.join('\n')
    """
    cat > config.yaml <<'YAMLEOF'
conditions:
${yaml_conditions}
samples:
${yaml_samples}
outdir: grouper_out
orphan: ${params.grouper_orphan ? 'True' : 'False'}
mincut: ${params.grouper_mincut ? 'True' : 'False'}
YAMLEOF

    Grouper --config config.yaml
    """
}
