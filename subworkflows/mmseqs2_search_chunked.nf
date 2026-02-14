/*
 * Subworkflow: MMSEQS2_SEARCH_CHUNKED
 *
 * Splits ORF queries into chunks, runs MMseqs2 searches in parallel, and merges results.
 * Provides same output as MMSEQS2_SEARCH but with 5-8× speedup through parallelization.
 */

include { SPLIT_ORFS          } from '../modules/mmseqs2_search_chunked'
include { MMSEQS2_SEARCH_CHUNK } from '../modules/mmseqs2_search_chunked'
include { MERGE_M8_RESULTS    } from '../modules/mmseqs2_search_chunked'

workflow MMSEQS2_SEARCH_CHUNKED {
    take:
    query_pep    // path: ORF peptide sequences
    db_path      // val: database path
    tag_name     // val: search tag (e.g., 'swissprot', 'pfam')

    main:
    // Step 1: Split ORFs into chunks
    SPLIT_ORFS(query_pep, tag_name)

    // Step 2: Create channel of chunks with indices
    ch_chunks = SPLIT_ORFS.out.chunks
        .flatten()
        .map { chunk_file ->
            def chunk_idx = (chunk_file.name =~ /orf_chunk_(\d+)\.pep/)[0][1]
            tuple(chunk_idx, chunk_file)
        }

    // Step 3: Run search on each chunk in parallel
    MMSEQS2_SEARCH_CHUNK(
        ch_chunks,
        db_path,
        tag_name
    )

    // Step 4: Merge all chunk results
    MERGE_M8_RESULTS(
        MMSEQS2_SEARCH_CHUNK.out.m8.collect(),
        tag_name
    )

    emit:
    m8 = MERGE_M8_RESULTS.out.m8
}
