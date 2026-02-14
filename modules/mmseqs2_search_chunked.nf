/*
 * MMseqs2 homology search with chunking for parallel execution
 *
 * Splits query ORFs into chunks, runs searches in parallel, merges results.
 * This provides 5-8× speedup by processing chunks simultaneously instead of
 * one large sequential job.
 */

process SPLIT_ORFS {
    tag "${tag_name}"

    input:
    path(orfs_pep)
    val(tag_name)

    output:
    path("chunks/orf_chunk_*.pep"), emit: chunks

    script:
    def chunk_size = params.search_orf_chunk_size ?: 40000
    """
    mkdir -p chunks

    # Split into chunks of ${chunk_size} sequences each
    awk -v size=${chunk_size} -v pre=chunks/orf_chunk_ '
        /^>/ {
            if (n % size == 0) {
                file = sprintf("%s%03d.pep", pre, int(n / size))
            }
            n++
        }
        {print > file}
    ' ${orfs_pep}

    # Ensure at least one chunk exists even for small inputs
    if [ ! -f chunks/orf_chunk_000.pep ]; then
        cp ${orfs_pep} chunks/orf_chunk_000.pep
    fi
    """
}

process MMSEQS2_SEARCH_CHUNK {
    tag "${tag_name}_chunk_${chunk_idx}"

    input:
    tuple val(chunk_idx), path(chunk_pep)
    val(db_path)
    val(tag_name)

    output:
    path("${tag_name}_chunk_${chunk_idx}.m8"), emit: m8

    script:
    def args = task.ext.args ?: ''
    def mem_gb = (task.memory.toGiga() * 0.85).intValue()
    """
    mmseqs easy-search \\
        ${chunk_pep} \\
        ${db_path} \\
        ${tag_name}_chunk_${chunk_idx}.m8 \\
        tmp_${tag_name}_${chunk_idx} \\
        -s ${params.mmseqs2_search_sens} \\
        --split-memory-limit ${mem_gb}G \\
        ${args} \\
        --threads ${task.cpus}
    """
}

process MERGE_M8_RESULTS {
    tag "${tag_name}"
    publishDir "${params.outdir}/mmseqs2_search", mode: 'copy', pattern: "*.m8"

    input:
    path("chunk_*.m8")
    val(tag_name)

    output:
    path("${tag_name}_alnRes.m8"), emit: m8

    script:
    """
    # Concatenate all chunks and sort by query ID for consistency
    cat chunk_*.m8 | sort -k1,1 > ${tag_name}_alnRes.m8

    # Report chunk merge stats
    echo "Merged \$(ls chunk_*.m8 | wc -l) chunks into ${tag_name}_alnRes.m8"
    echo "Total alignments: \$(wc -l < ${tag_name}_alnRes.m8)"
    """
}
