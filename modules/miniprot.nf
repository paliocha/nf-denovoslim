/*
 * MINIPROT_ALIGN — protein-to-genome spliced alignment using miniprot.
 *
 * Aligns predicted proteins to a reference genome to validate gene
 * predictions.  The reference genome is passed as val (not path) to
 * avoid staging multi-GB files into the work directory — container
 * bind mounts make it accessible directly.
 *
 * Only runs when --reference_genome is provided.
 */

process MINIPROT_ALIGN {
    tag "${species_label}"

    input:
    path(proteins)
    val(reference_genome)
    val(species_label)

    output:
    path("miniprot_alignments.gff"), emit: gff

    script:
    """
    miniprot \\
        --gff \\
        -t ${task.cpus} \\
        "${reference_genome}" \\
        ${proteins} \\
        > miniprot_alignments.gff

    N_PROT=\$(grep -c '^>' ${proteins})
    N_ALN=\$(awk -F'\\t' '\$3 == "mRNA"' miniprot_alignments.gff | wc -l)
    N_MAPPED=\$(awk -F'\\t' '\$3 == "mRNA" {
        n = split(\$9, attrs, ";")
        for (i = 1; i <= n; i++) {
            if (attrs[i] ~ /^Target=/) {
                split(attrs[i], parts, "=")
                split(parts[2], tparts, " ")
                print tparts[1]
            }
        }
    }' miniprot_alignments.gff | sort -u | wc -l)
    echo "miniprot: \$N_PROT proteins, \$N_ALN alignments, \$N_MAPPED mapped"
    """
}

/*
 * MINIPROT_FILTER — Filter proteins based on genome validation.
 *
 * Keeps proteins that either:
 *   1. Map to the genome with Identity >= min_identity, OR
 *   2. Are >= min_length amino acids (trusted by size alone)
 *
 * Short proteins that lack genomic support are removed as potential
 * prediction artifacts.  Only runs when --reference_genome is provided.
 */

process MINIPROT_FILTER {
    tag "${species_label}"

    input:
    path(proteins)
    path(gff)
    val(species_label)

    output:
    path("${species_label}.faa"), emit: faa
    path("miniprot_stats.txt"),   emit: stats

    script:
    def min_identity = params.miniprot_min_identity ?: 0.5
    def min_length   = params.miniprot_min_length   ?: 75
    """
    miniprot_filter.py \\
        --proteins ${proteins} \\
        --gff ${gff} \\
        --min-identity ${min_identity} \\
        --min-length ${min_length} \\
        --out ${species_label}.faa \\
        --stats miniprot_stats.txt
    """
}
