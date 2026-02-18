/*
 * Subworkflow: FRAMESHIFT_CORRECTION_CHUNKED
 *
 * Splits SuperTranscripts into chunks, runs Diamond blastx + frameshift
 * correction in parallel, and merges results. Provides same output as
 * FRAMESHIFT_CORRECTION but with 3-5x speedup (8h -> ~2h).
 *
 * Diamond and Python correction are separate processes so that Diamond
 * results stay cached if the Python step needs to be re-run.
 */

include { SPLIT_SUPERTRANSCRIPTS_FS    } from '../modules/frameshift_correction_chunked'
include { DIAMOND_BLASTX_CHUNK         } from '../modules/frameshift_correction_chunked'
include { CORRECT_FRAMESHIFTS_CHUNK    } from '../modules/frameshift_correction_chunked'
include { MERGE_FRAMESHIFT_RESULTS     } from '../modules/frameshift_correction_chunked'

workflow FRAMESHIFT_CORRECTION_CHUNKED {
    take:
    supertranscripts_fasta    // path: input FASTA file
    diamond_db                // val: path to Diamond DB
    species_label             // val: species label for tags

    main:
    // Step 1: Split SuperTranscripts into chunks
    SPLIT_SUPERTRANSCRIPTS_FS(supertranscripts_fasta, species_label)

    // Step 2: Create channel of chunks with indices
    ch_chunks = SPLIT_SUPERTRANSCRIPTS_FS.out.chunks
        .flatten()
        .map { chunk_file ->
            def chunk_idx = (chunk_file.name =~ /fs_chunk_(\d+)\.fasta/)[0][1]
            tuple(chunk_idx, chunk_file)
        }

    // Step 3: Run Diamond blastx on each chunk (heavy computation, cached independently)
    DIAMOND_BLASTX_CHUNK(ch_chunks, diamond_db, species_label)

    // Step 4: Apply frameshift corrections using Python (lightweight, separate cache)
    CORRECT_FRAMESHIFTS_CHUNK(DIAMOND_BLASTX_CHUNK.out.results, species_label)

    // Step 5: Merge all chunk results
    MERGE_FRAMESHIFT_RESULTS(
        CORRECT_FRAMESHIFTS_CHUNK.out.fasta.map { idx, fasta -> fasta }.collect(),
        CORRECT_FRAMESHIFTS_CHUNK.out.stats.collect(),
        species_label
    )

    emit:
    fasta = MERGE_FRAMESHIFT_RESULTS.out.fasta
    stats = MERGE_FRAMESHIFT_RESULTS.out.stats
}
