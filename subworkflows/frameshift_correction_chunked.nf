/*
 * Subworkflow: FRAMESHIFT_CORRECTION_CHUNKED
 *
 * Splits SuperTranscripts into chunks, runs Diamond blastx + frameshift
 * correction in parallel, and merges results. Provides same output as
 * FRAMESHIFT_CORRECTION but with 3-5× speedup (8h → ~2h).
 */

include { SPLIT_SUPERTRANSCRIPTS_FS    } from '../modules/frameshift_correction_chunked'
include { FRAMESHIFT_CORRECTION_CHUNK  } from '../modules/frameshift_correction_chunked'
include { MERGE_FRAMESHIFT_RESULTS     } from '../modules/frameshift_correction_chunked'

workflow FRAMESHIFT_CORRECTION_CHUNKED {
    take:
    supertranscripts_fasta    // path: input FASTA file

    main:
    // Step 1: Split SuperTranscripts into chunks
    SPLIT_SUPERTRANSCRIPTS_FS(supertranscripts_fasta)

    // Step 2: Create channel of chunks with indices
    ch_chunks = SPLIT_SUPERTRANSCRIPTS_FS.out.chunks
        .flatten()
        .map { chunk_file ->
            def chunk_idx = (chunk_file.name =~ /fs_chunk_(\d+)\.fasta/)[0][1]
            tuple(chunk_idx, chunk_file)
        }

    // Step 3: Run Diamond + frameshift correction on each chunk in parallel
    FRAMESHIFT_CORRECTION_CHUNK(ch_chunks)

    // Step 4: Merge all chunk results
    MERGE_FRAMESHIFT_RESULTS(
        FRAMESHIFT_CORRECTION_CHUNK.out.fasta.map { idx, fasta -> fasta }.collect(),
        FRAMESHIFT_CORRECTION_CHUNK.out.stats.collect()
    )

    emit:
    fasta = MERGE_FRAMESHIFT_RESULTS.out.fasta
    stats = MERGE_FRAMESHIFT_RESULTS.out.stats
}
