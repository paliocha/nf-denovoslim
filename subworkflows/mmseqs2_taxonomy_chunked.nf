/*
 * Subworkflow: MMSEQS2_TAXONOMY_CHUNKED
 *
 * Splits SuperTranscripts into chunks, runs taxonomy classification and filtering
 * in parallel, and merges results. Provides same output as MMSEQS2_TAXONOMY but
 * with 5-10× speedup through parallelization (24h → 3-4h).
 */

include { SPLIT_SUPERTRANSCRIPTS   } from '../modules/mmseqs2_taxonomy_chunked'
include { MMSEQS2_TAXONOMY_CHUNK   } from '../modules/mmseqs2_taxonomy_chunked'
include { MERGE_TAXONOMY_RESULTS   } from '../modules/mmseqs2_taxonomy_chunked'

workflow MMSEQS2_TAXONOMY_CHUNKED {
    take:
    supertranscripts_fasta    // path: input FASTA file

    main:
    // Step 1: Split SuperTranscripts into chunks
    SPLIT_SUPERTRANSCRIPTS(supertranscripts_fasta)

    // Step 2: Create channel of chunks with indices
    ch_chunks = SPLIT_SUPERTRANSCRIPTS.out.chunks
        .flatten()
        .map { chunk_file ->
            def chunk_idx = (chunk_file.name =~ /st_chunk_(\d+)\.fasta/)[0][1]
            tuple(chunk_idx, chunk_file)
        }

    // Step 3: Run taxonomy classification + filtering on each chunk in parallel
    MMSEQS2_TAXONOMY_CHUNK(ch_chunks)

    // Step 4: Merge all chunk results
    MERGE_TAXONOMY_RESULTS(
        MMSEQS2_TAXONOMY_CHUNK.out.fasta.map { idx, fasta -> fasta }.collect(),
        MMSEQS2_TAXONOMY_CHUNK.out.tsv.map { idx, tsv -> tsv }.collect(),
        MMSEQS2_TAXONOMY_CHUNK.out.stats.collect(),
        supertranscripts_fasta
    )

    emit:
    fasta   = MERGE_TAXONOMY_RESULTS.out.fasta
    lca_tsv = MERGE_TAXONOMY_RESULTS.out.lca_tsv
    stats   = MERGE_TAXONOMY_RESULTS.out.stats
}
