/*
 * SortMeRNA — rRNA filtering (matching nf-core/rnaseq v3.22.2)
 */

process SORTMERNA_INDEX {
    label 'process_medium'
    tag 'sortmerna_index'

    input:
    path(fastas)

    output:
    path("sortmerna_idx"), emit: index

    script:
    def refs = fastas.collect { "--ref ${it}" }.join(' ')
    """
    sortmerna \\
        ${refs} \\
        --workdir . \\
        --index 1

    mkdir -p sortmerna_idx
    mv idx/ sortmerna_idx/ 2>/dev/null || true
    mv kvdb/ sortmerna_idx/ 2>/dev/null || true
    """
}

process SORTMERNA {
    label 'process_high'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)
    path(fastas)
    path(index)

    output:
    tuple val(sample_id), path("non_rRNA_reads_fwd.fq.gz"), path("non_rRNA_reads_rev.fq.gz"), emit: reads
    tuple val(sample_id), path("sortmerna.log"), emit: log

    script:
    def refs = fastas.collect { "--ref ${it}" }.join(' ')
    """
    # Symlink pre-built index so SortMeRNA finds it at ./idx/
    ln -s ${index}/idx idx

    sortmerna \\
        ${refs} \\
        --reads ${reads_1} --reads ${reads_2} \\
        --threads ${task.cpus} \\
        --workdir . \\
        --aligned rRNA_reads \\
        --fastx \\
        --other non_rRNA_reads \\
        --paired_in \\
        --out2 \\
        --num_alignments 1 \\
        -v \\
        --index 0 \\
        2>&1 | tee sortmerna.log

    # SortMeRNA outputs non_rRNA_reads_fwd.f[qa][.gz] and _rev
    # Ensure gzipped output
    if [ -f non_rRNA_reads_fwd.fq ] && [ ! -f non_rRNA_reads_fwd.fq.gz ]; then
        gzip non_rRNA_reads_fwd.fq
        gzip non_rRNA_reads_rev.fq
    fi
    """
}
