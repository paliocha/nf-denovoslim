#!/usr/bin/env nextflow

/*
 * nf-denovoslim — collapse a Trinity de novo transcriptome into a
 * non-redundant gene set with SuperTranscripts, one best protein per
 * gene, gene-level Salmon quantification, and functional annotation.
 */

// --- Module imports ---

include { SORTMERNA_INDEX                  } from './modules/sortmerna'
include { SORTMERNA                        } from './modules/sortmerna'
include { SALMON_INDEX as SALMON_INDEX_INITIAL } from './modules/salmon_index'
include { SALMON_INDEX as SALMON_INDEX_FINAL   } from './modules/salmon_index'
include { SALMON_QUANT as SALMON_QUANT_INITIAL } from './modules/salmon_quant'
include { SALMON_QUANT as SALMON_QUANT_FINAL   } from './modules/salmon_quant'
include { CORSET                           } from './modules/corset'
include { LACE                             } from './modules/lace'
include { MMSEQS2_TAXONOMY                 } from './modules/mmseqs2_taxonomy'
include { DIAMOND_BLASTX                   } from './modules/frameshift_correction'
include { CORRECT_FRAMESHIFTS              } from './modules/frameshift_correction'
include { TD2_LONGORFS                     } from './modules/td2_longorfs'
include { MMSEQS2_SEARCH as MMSEQS2_SEARCH_SWISSPROT } from './modules/mmseqs2_search'
include { MMSEQS2_SEARCH as MMSEQS2_SEARCH_PFAM      } from './modules/mmseqs2_search'
include { TD2_PREDICT                      } from './modules/td2_predict'
include { SELECT_BEST_ORF                  } from './modules/select_best_orf'
include { VALIDATE_IDS                     } from './modules/validate_ids'
include { BUSCO as BUSCO_TRINITY           } from './modules/busco'
include { BUSCO as BUSCO_QC                } from './modules/busco'
include { TRANSANNOT                       } from './modules/transannot'
include { THINNING_REPORT                  } from './modules/thinning_report'

// --- Typed pipeline parameters (Nextflow ≥25.10, strict syntax) ---

params {
    // --- Input (required — no default) ---
    trinity_fasta:       Path
    samplesheet:         Path
    species_label:       String  = 'species_X'

    // --- SortMeRNA rRNA databases ---
    sortmerna_db_urls:   List    = [
        'https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/rfam-5.8s-database-id98.fasta',
        'https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/rfam-5s-database-id98.fasta',
        'https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-arc-16s-id95.fasta',
        'https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-arc-23s-id98.fasta',
        'https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-bac-16s-id90.fasta',
        'https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-bac-23s-id98.fasta',
        'https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-euk-18s-id95.fasta',
        'https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-euk-28s-id98.fasta'
    ]
    sortmerna_db_dir:    String? = null   // Pre-downloaded dir; if null, downloads from URLs

    // --- Databases (pre-built MMseqs2 — set via CLI or site config) ---
    mmseqs2_swissprot:   Path
    mmseqs2_pfam:        Path
    mmseqs2_eggnog:      Path
    mmseqs2_taxonomy_db: Path
    eggnog_annotations:  Path
    busco_lineage:       String              // e.g. 'poales_odb12', 'eudicots_odb12'
    diamond_db:          Path

    // --- Cluster (optional) ---
    unix_group:          String? = null
    orion_exclude_nodes: String? = null      // e.g. 'cn-37' (Orion-specific)

    // --- Search params ---
    mmseqs2_search_sens: Float   = 7.0

    // --- Taxonomy filter (NCBI taxon ID — includes all descendants) ---
    filter_taxon:        Integer = 33090     // Viridiplantae (green plants)

    // --- TD2 params ---
    td2_min_orf_length:  Integer = 90
    td2_strand_specific: Boolean = true

    // --- Output ---
    outdir:              String  = './results'
}

// --- Main workflow ---

workflow {

    main:

    // Parse samplesheet: sample,fastq_1,fastq_2,strandedness[,condition]
    ch_samplesheet = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def condition = row.condition ?: Utils.extractCondition(sample_id)
            def reads_1   = file(row.fastq_1, checkIfExists: true)
            def reads_2   = file(row.fastq_2, checkIfExists: true)
            [ sample_id, condition, reads_1, reads_2 ]
        }

    ch_trinity = channel.fromPath(params.trinity_fasta, checkIfExists: true)

    // BUSCO baseline on raw Trinity (transcriptome mode, no dependencies)
    BUSCO_TRINITY(ch_trinity, params.species_label, 'trinity')

    // -- SortMeRNA: rRNA filtering --

    if (params.sortmerna_db_dir) {
        ch_sortmerna_fastas = channel
            .fromPath("${params.sortmerna_db_dir}/*.fasta")
            .collect()
    } else {
        ch_sortmerna_fastas = channel
            .fromList(params.sortmerna_db_urls)
            .map { url -> file(url) }
            .collect()
    }

    SORTMERNA_INDEX(ch_sortmerna_fastas)

    ch_reads_for_sortmerna = ch_samplesheet
        .map { sample_id, _condition, r1, r2 -> [ sample_id, r1, r2 ] }

    SORTMERNA(
        ch_reads_for_sortmerna,
        ch_sortmerna_fastas,
        SORTMERNA_INDEX.out.index
    )

    // Rejoin filtered reads with condition metadata (re-read samplesheet
    // to avoid consuming the queue channel used for reads)
    ch_condition_map = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> [ row.sample, row.condition ?: Utils.extractCondition(row.sample) ] }
        .toList()
        .map { list -> list.collectEntries { entry -> [ entry[0], entry[1] ] } }

    ch_filtered_reads = SORTMERNA.out.reads
        .combine(ch_condition_map)
        .map { sample_id, r1, r2, cmap ->
            [ sample_id, cmap[sample_id] ?: Utils.extractCondition(sample_id), r1, r2 ]
        }

    // -- Salmon initial quant (full Trinity, --dumpEq for Corset) --

    SALMON_INDEX_INITIAL(ch_trinity)

    ch_reads_for_salmon = ch_filtered_reads
        .map { sample_id, _condition, r1, r2 -> [ sample_id, r1, r2 ] }

    SALMON_QUANT_INITIAL(
        ch_reads_for_salmon,
        SALMON_INDEX_INITIAL.out.index.first(),
        'quant'
    )

    // -- Corset: transcript-to-gene clustering --

    ch_all_quants = SALMON_QUANT_INITIAL.out.quant_dir.collect()

    // Build Corset -g/-n from a fresh samplesheet read (avoids queue-channel fork)
    ch_sample_conditions = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            [ sample_id: row.sample, condition: row.condition ?: Utils.extractCondition(row.sample) ]
        }
        .toSortedList { a, b ->
            a.condition <=> b.condition ?: a.sample_id <=> b.sample_id
        }

    CORSET(ch_all_quants, ch_sample_conditions, params.species_label)

    // -- Lace: build SuperTranscripts --

    LACE(ch_trinity, CORSET.out.clust, params.species_label)

    // -- Taxonomy filter (keep Viridiplantae) --

    MMSEQS2_TAXONOMY(
        LACE.out.fasta,
        params.mmseqs2_taxonomy_db,
        params.mmseqs2_search_sens,
        params.filter_taxon,
        params.species_label
    )

    // -- Frameshift correction --

    DIAMOND_BLASTX(MMSEQS2_TAXONOMY.out.fasta, params.diamond_db, params.species_label)
    CORRECT_FRAMESHIFTS(MMSEQS2_TAXONOMY.out.fasta, DIAMOND_BLASTX.out.tsv, params.species_label)

    // -- ORF prediction (TD2 + homology support) --

    TD2_LONGORFS(CORRECT_FRAMESHIFTS.out.fasta, params.species_label)

    MMSEQS2_SEARCH_SWISSPROT(
        TD2_LONGORFS.out.longest_orfs_pep,
        params.mmseqs2_swissprot,
        'swissprot'
    )

    MMSEQS2_SEARCH_PFAM(
        TD2_LONGORFS.out.longest_orfs_pep,
        params.mmseqs2_pfam,
        'pfam'
    )

    TD2_PREDICT(
        CORRECT_FRAMESHIFTS.out.fasta,
        MMSEQS2_SEARCH_SWISSPROT.out.m8,
        MMSEQS2_SEARCH_PFAM.out.m8,
        TD2_LONGORFS.out.td2_dir,
        params.species_label
    )

    // -- Best ORF selection --

    SELECT_BEST_ORF(
        TD2_PREDICT.out.psauron_scores,
        TD2_PREDICT.out.pep,
        TD2_PREDICT.out.gff3,
        params.species_label
    )

    // -- Gene-level Salmon quant on SuperTranscripts --

    SALMON_INDEX_FINAL(CORRECT_FRAMESHIFTS.out.fasta)

    SALMON_QUANT_FINAL(
        ch_reads_for_salmon,
        SALMON_INDEX_FINAL.out.index.first(),
        'st_quant'
    )

    // -- ID validation --

    ch_first_quant_sf = SALMON_QUANT_FINAL.out.quant_dir
        .first()
        .map { qdir -> file("${qdir}/quant.sf") }

    VALIDATE_IDS(
        ch_first_quant_sf,
        SELECT_BEST_ORF.out.faa,
        params.species_label
    )

    // -- BUSCO QC + TransAnnot (parallel) --

    BUSCO_QC(SELECT_BEST_ORF.out.faa, params.species_label, 'final')

    TRANSANNOT(
        SELECT_BEST_ORF.out.faa,
        file(params.eggnog_annotations, checkIfExists: true),
        params.species_label,
        params.mmseqs2_pfam,
        params.mmseqs2_eggnog,
        params.mmseqs2_swissprot
    )

    // -- Thinning report --

    THINNING_REPORT(
        ch_trinity,
        CORRECT_FRAMESHIFTS.out.fasta,
        CORSET.out.clust,
        SELECT_BEST_ORF.out.map,
        SELECT_BEST_ORF.out.faa,
        SALMON_QUANT_INITIAL.out.quant_dir.collect(),
        SALMON_QUANT_FINAL.out.quant_dir.collect(),
        BUSCO_TRINITY.out.outdir,
        BUSCO_QC.out.outdir,
        VALIDATE_IDS.out.report,
        SORTMERNA.out.log.map { _sample_id, logfile -> logfile }.collect(),
        MMSEQS2_TAXONOMY.out.breakdown,
        params.species_label
    )

    // -- Published outputs (structure matches README § Output) --

    publish:
    clustering_clusters = CORSET.out.clust
    clustering_counts   = CORSET.out.counts
    supertranscripts    = LACE.out.fasta
    taxonomy_fasta      = MMSEQS2_TAXONOMY.out.fasta
    taxonomy_lca        = MMSEQS2_TAXONOMY.out.lca_tsv
    taxonomy_stats      = MMSEQS2_TAXONOMY.out.stats
    taxonomy_breakdown  = MMSEQS2_TAXONOMY.out.breakdown
    frameshift_stats    = CORRECT_FRAMESHIFTS.out.stats
    search_swissprot    = MMSEQS2_SEARCH_SWISSPROT.out.m8
    search_pfam         = MMSEQS2_SEARCH_PFAM.out.m8
    proteome            = SELECT_BEST_ORF.out.faa
    best_orf_gff3       = SELECT_BEST_ORF.out.gff3
    orf_gene_map        = SELECT_BEST_ORF.out.map
    salmon_quant        = SALMON_QUANT_FINAL.out.quant_dir
    busco_trinity       = BUSCO_TRINITY.out.outdir
    busco_qc            = BUSCO_QC.out.outdir
    transannot          = TRANSANNOT.out.annotation
    thinning_report     = THINNING_REPORT.out.report

    onComplete:
    log.info ""
    log.info "Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration          : ${workflow.duration}"
    log.info "Output dir        : ${params.outdir}"
    log.info ""
}

// -- Workflow output definitions (Nextflow ≥25.10) --
// Maps published channels to output directory structure.
// Default outputDir is params.outdir (see nextflow.config).

output {
    clustering_clusters { path 'clustering' }
    clustering_counts   { path 'clustering' }
    supertranscripts    { path 'supertranscripts' }
    taxonomy_fasta      { path 'taxonomy' }
    taxonomy_lca        { path 'taxonomy' }
    taxonomy_stats      { path 'taxonomy' }
    taxonomy_breakdown  { path 'taxonomy' }
    frameshift_stats    { path 'frameshift_correction' }
    search_swissprot    { path 'mmseqs2_search' }
    search_pfam         { path 'mmseqs2_search' }
    proteome            { path 'proteins' }
    best_orf_gff3       { path 'annotation' }
    orf_gene_map        { path 'annotation' }
    salmon_quant        { path 'salmon_final' }
    busco_trinity       { path 'qc' }
    busco_qc            { path 'qc' }
    transannot          { path 'transannot' }
    thinning_report     { path '.' }
}
