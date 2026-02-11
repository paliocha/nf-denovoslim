/*
 * SuperTranscripts — merge transcripts per gene group
 */

process SUPERTRANSCRIPTS {
    label 'process_high'
    tag "${params.species_label}"

    publishDir "${params.outdir}/supertranscripts", mode: 'copy'

    input:
    path(deduped_fasta)
    path(grouper_clust)

    output:
    path("trinity_genes.fasta"), emit: fasta

    script:
    """
    # Convert Grouper mag.flat.clust to Trinity gene_trans_map format
    # mag.flat.clust: transcript_id\\tcluster_id
    # gene_trans_map: gene_id\\ttranscript_id  (columns swapped)
    awk '{print \$2"\\t"\$1}' ${grouper_clust} > gene_trans_map.tsv

    python \$TRINITY_HOME/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \\
        --trinity_fasta ${deduped_fasta} \\
        --gene_trans_map gene_trans_map.tsv

    # The script outputs trinity_genes.fasta in the working directory
    """
}
