#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-denovoslim — Trinity Assembly Thinning Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Collapses a fragmented Trinity de novo transcriptome assembly into a non-redundant
    gene set with SuperTranscripts, gene-level Salmon quantification, and one best
    protein per gene.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ──────────────────────────────────────────────────────────────────────────────
//  Module imports (DSL2 aliasing for reusable processes)
// ──────────────────────────────────────────────────────────────────────────────

include { SORTMERNA_INDEX                            } from './modules/sortmerna'
include { SORTMERNA                                  } from './modules/sortmerna'
include { MMSEQS2_CLUSTER_NT                         } from './modules/mmseqs2_cluster_nt'
include { SALMON_INDEX as SALMON_INDEX_INITIAL       } from './modules/salmon_index'
include { SALMON_INDEX as SALMON_INDEX_FINAL         } from './modules/salmon_index'
include { SALMON_QUANT as SALMON_QUANT_INITIAL       } from './modules/salmon_quant'
include { SALMON_QUANT as SALMON_QUANT_FINAL         } from './modules/salmon_quant'
include { CORSET                                     } from './modules/corset'
include { LACE                                       } from './modules/lace'
include { MMSEQS2_TAXONOMY                            } from './modules/mmseqs2_taxonomy'
include { DIAMOND_BLASTX                              } from './modules/frameshift_correction'
include { CORRECT_FRAMESHIFTS                          } from './modules/frameshift_correction'
include { TD2_LONGORFS                                               } from './modules/td2_longorfs'
include { MMSEQS2_SEARCH_CHUNKED as MMSEQS2_SEARCH_CHUNKED_SWISSPROT } from './subworkflows/mmseqs2_search_chunked'
include { MMSEQS2_SEARCH_CHUNKED as MMSEQS2_SEARCH_CHUNKED_PFAM      } from './subworkflows/mmseqs2_search_chunked'
include { TD2_PREDICT                                                } from './modules/td2_predict'
include { SELECT_BEST_ORF                            } from './modules/select_best_orf'
include { VALIDATE_IDS                               } from './modules/validate_ids'
include { BUSCO_QC                                   } from './modules/busco'
include { TRANSANNOT                                 } from './modules/transannot'
include { THINNING_REPORT                            } from './modules/thinning_report'

// ──────────────────────────────────────────────────────────────────────────────
//  Input validation
// ──────────────────────────────────────────────────────────────────────────────

if (!params.trinity_fasta)       { error "Please provide --trinity_fasta" }
if (!params.samplesheet)         { error "Please provide --samplesheet" }
if (!params.mmseqs2_swissprot)   { error "Please provide --mmseqs2_swissprot (path to MMseqs2 SwissProt DB)" }
if (!params.mmseqs2_pfam)        { error "Please provide --mmseqs2_pfam (path to MMseqs2 Pfam DB)" }
if (!params.mmseqs2_eggnog)      { error "Please provide --mmseqs2_eggnog (path to MMseqs2 eggNOG DB)" }
if (!params.mmseqs2_taxonomy_db) { error "Please provide --mmseqs2_taxonomy_db (path to MMseqs2 taxonomy DB)" }
if (!params.eggnog_annotations)  { error "Please provide --eggnog_annotations (path to eggNOG annotation TSV)" }
if (!params.busco_lineage)       { error "Please provide --busco_lineage (e.g. 'poales_odb12', 'eudicots_odb12')" }

// ──────────────────────────────────────────────────────────────────────────────
//  Helper: extract condition from sample name (fallback if no 'condition' column)
//  e.g. BMAX56_T4_L -> T4_L
// ──────────────────────────────────────────────────────────────────────────────

def extractCondition(sample_name) {
    def parts = sample_name.split('_')
    if (parts.size() < 3) {
        log.warn "Cannot extract condition from sample '${sample_name}' — using full name as condition"
        return sample_name
    }
    // Last two parts are Timepoint and Tissue
    return "${parts[-2]}_${parts[-1]}"
}

// ──────────────────────────────────────────────────────────────────────────────
//  Main workflow
// ──────────────────────────────────────────────────────────────────────────────

workflow {

    // --- Parse samplesheet (nf-core/rnaseq format) ---
    // CSV: sample,fastq_1,fastq_2,strandedness[,condition]
    // If 'condition' column is present, it is used for Corset grouping;
    // otherwise the condition is extracted from sample name (last two _-separated parts).
    ch_samplesheet = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def condition = row.condition ?: extractCondition(sample_id)
            def reads_1   = file(row.fastq_1, checkIfExists: true)
            def reads_2   = file(row.fastq_2, checkIfExists: true)
            [ sample_id, condition, reads_1, reads_2 ]
        }

    // --- Trinity assembly ---
    ch_trinity = Channel.fromPath(params.trinity_fasta, checkIfExists: true)

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 0: SortMeRNA — rRNA filtering                               ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    // Prepare rRNA databases
    if (params.sortmerna_db_dir) {
        ch_sortmerna_fastas = Channel
            .fromPath("${params.sortmerna_db_dir}/*.fasta")
            .collect()
    } else {
        ch_sortmerna_fastas = Channel
            .fromList(params.sortmerna_db_urls)
            .map { url -> file(url) }
            .collect()
    }

    // Build index once
    SORTMERNA_INDEX(ch_sortmerna_fastas)

    // Filter each sample
    ch_reads_for_sortmerna = ch_samplesheet
        .map { sample_id, condition, r1, r2 -> [ sample_id, r1, r2 ] }

    SORTMERNA(
        ch_reads_for_sortmerna,
        ch_sortmerna_fastas,
        SORTMERNA_INDEX.out.index
    )

    // Recombine filtered reads with condition metadata
    // Derive condition from samplesheet (with condition column support)
    // Re-read samplesheet to build a lookup map for condition per sample
    ch_condition_map = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> [ row.sample, row.condition ?: extractCondition(row.sample) ] }
        .toList()
        .map { list -> list.collectEntries { [ it[0], it[1] ] } }

    ch_filtered_reads = SORTMERNA.out.reads
        .combine(ch_condition_map)
        .map { sample_id, r1, r2, cmap ->
            [ sample_id, cmap[sample_id] ?: extractCondition(sample_id), r1, r2 ]
        }

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 1: MMseqs2 nucleotide clustering (97% identity dedup)        ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    MMSEQS2_CLUSTER_NT(ch_trinity)

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 2-3: Salmon initial index + quant (with --dumpEq for Corset) ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    SALMON_INDEX_INITIAL(MMSEQS2_CLUSTER_NT.out.rep_fasta)

    ch_reads_for_salmon = ch_filtered_reads
        .map { sample_id, condition, r1, r2 -> [ sample_id, r1, r2 ] }

    SALMON_QUANT_INITIAL(
        ch_reads_for_salmon,
        SALMON_INDEX_INITIAL.out.index.first(),
        'quant'
    )

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 4: Corset — hierarchical transcript-to-gene clustering       ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    // Collect all quant directories
    ch_all_quants = SALMON_QUANT_INITIAL.out.quant_dir.collect()

    // Build sample-condition metadata for Corset -g/-n flags
    // Derived from samplesheet directly (not ch_filtered_reads) to avoid
    // DSL2 queue-channel fork that splits items between consumers
    ch_sample_conditions = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            [ sample_id: row.sample, condition: row.condition ?: extractCondition(row.sample) ]
        }
        .toSortedList { a, b ->
            a.condition <=> b.condition ?: a.sample_id <=> b.sample_id
        }

    CORSET(ch_all_quants, ch_sample_conditions, params.species_label)

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 5: Lace — build SuperTranscripts from Corset clusters        ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    LACE(
        MMSEQS2_CLUSTER_NT.out.rep_fasta,
        CORSET.out.clust,
        params.species_label
    )

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 5b: Taxonomy filter — keep only Streptophyta SuperTranscripts║
    // ║  (Single process on node-local SSD for fast DB I/O)                ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    MMSEQS2_TAXONOMY(
        LACE.out.fasta,
        params.mmseqs2_taxonomy_db,
        params.mmseqs2_search_sens,
        params.filter_taxon,
        params.species_label
    )

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 5c: Frameshift correction — fix assembly frameshifts         ║
    // ║  (Diamond DB copied to node-local SSD; correction cached separately)║
    // ╚══════════════════════════════════════════════════════════════════════╝

    DIAMOND_BLASTX(MMSEQS2_TAXONOMY.out.fasta, params.diamond_db, params.species_label)
    CORRECT_FRAMESHIFTS(MMSEQS2_TAXONOMY.out.fasta, DIAMOND_BLASTX.out.tsv, params.species_label)

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 6-9: TD2 ORF prediction with homology support                ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    TD2_LONGORFS(CORRECT_FRAMESHIFTS.out.fasta, params.species_label)

    // Steps 7 & 8 run in parallel with chunking for 5-8× speedup
    // (DB paths passed as val — no staging of multi-GB DBs)
    MMSEQS2_SEARCH_CHUNKED_SWISSPROT(
        TD2_LONGORFS.out.longest_orfs_pep,
        params.mmseqs2_swissprot,
        'swissprot'
    )

    MMSEQS2_SEARCH_CHUNKED_PFAM(
        TD2_LONGORFS.out.longest_orfs_pep,
        params.mmseqs2_pfam,
        'pfam'
    )

    // Step 9: TD2.Predict with combined homology hits
    TD2_PREDICT(
        CORRECT_FRAMESHIFTS.out.fasta,
        MMSEQS2_SEARCH_CHUNKED_SWISSPROT.out.m8,
        MMSEQS2_SEARCH_CHUNKED_PFAM.out.m8,
        TD2_LONGORFS.out.td2_dir,
        params.species_label
    )

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 10: Select best ORF per gene                                 ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    SELECT_BEST_ORF(
        TD2_PREDICT.out.psauron_scores,
        TD2_PREDICT.out.pep,
        TD2_PREDICT.out.gff3,
        params.species_label
    )

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 11: Salmon final quant on SuperTranscripts (gene-level)      ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    SALMON_INDEX_FINAL(CORRECT_FRAMESHIFTS.out.fasta)

    SALMON_QUANT_FINAL(
        ch_reads_for_salmon,
        SALMON_INDEX_FINAL.out.index.first(),
        'st_quant'
    )

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 12: Validate ID consistency                                  ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    // Use the first sample's quant.sf for ID validation
    ch_first_quant_sf = SALMON_QUANT_FINAL.out.quant_dir
        .first()
        .map { qdir -> file("${qdir}/quant.sf") }

    VALIDATE_IDS(
        ch_first_quant_sf,
        SELECT_BEST_ORF.out.faa,
        params.species_label
    )

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 13-14: BUSCO QC + TransAnnot (run in parallel)               ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    BUSCO_QC(SELECT_BEST_ORF.out.faa, params.species_label)

    TRANSANNOT(
        SELECT_BEST_ORF.out.faa,
        file(params.eggnog_annotations, checkIfExists: true),
        params.species_label,
        params.mmseqs2_pfam,
        params.mmseqs2_eggnog,
        params.mmseqs2_swissprot
    )

    // ╔══════════════════════════════════════════════════════════════════════╗
    // ║  STEP 15: Thinning report                                          ║
    // ╚══════════════════════════════════════════════════════════════════════╝

    THINNING_REPORT(
        ch_trinity,
        MMSEQS2_CLUSTER_NT.out.rep_fasta,
        CORRECT_FRAMESHIFTS.out.fasta,
        CORSET.out.clust,
        SELECT_BEST_ORF.out.map,
        SELECT_BEST_ORF.out.faa,
        SALMON_QUANT_INITIAL.out.quant_dir.collect(),
        SALMON_QUANT_FINAL.out.quant_dir.collect(),
        BUSCO_QC.out.outdir,
        VALIDATE_IDS.out.report,
        SORTMERNA.out.log.map { sample_id, log -> log }.collect(),
        params.species_label
    )
}

// ──────────────────────────────────────────────────────────────────────────────
//  Completion handler
// ──────────────────────────────────────────────────────────────────────────────

workflow.onComplete {
    log.info ""
    log.info "Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration          : ${workflow.duration}"
    log.info "Output dir        : ${params.outdir}"
    log.info ""
}
