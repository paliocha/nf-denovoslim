#!/usr/bin/env nextflow

/*
 * nf-denovoslim — collapse a Trinity de novo transcriptome into a
 * non-redundant gene set with representative transcripts, merged
 * TD2+MetaEuk+GeneMarkS-T proteins, gene-level Salmon quantification,
 * and functional annotation.
 */

// --- Module imports ---

include { SORTMERNA_INDEX                  } from './modules/sortmerna'
include { SORTMERNA                        } from './modules/sortmerna'
include { SALMON_INDEX as SALMON_INDEX_INITIAL } from './modules/salmon_index'
include { SALMON_INDEX as SALMON_INDEX_FINAL   } from './modules/salmon_index'
include { SALMON_QUANT as SALMON_QUANT_INITIAL } from './modules/salmon_quant'
include { SALMON_QUANT as SALMON_QUANT_FINAL   } from './modules/salmon_quant'
include { CORSET                           } from './modules/corset'
include { SELECT_REP                       } from './modules/select_representative'
include { MMSEQS2_CLUSTER                  } from './modules/mmseqs2_cluster'
include { MMSEQS2_TAXONOMY                 } from './modules/mmseqs2_taxonomy'
include { DIAMOND_BLASTX                   } from './modules/frameshift_correction'
include { CORRECT_FRAMESHIFTS              } from './modules/frameshift_correction'
include { TD2_LONGORFS                     } from './modules/td2_longorfs'
include { MMSEQS2_SEARCH as MMSEQS2_SEARCH_SWISSPROT } from './modules/mmseqs2_search'
include { MMSEQS2_SEARCH as MMSEQS2_SEARCH_PFAM      } from './modules/mmseqs2_search'
include { TD2_PREDICT                      } from './modules/td2_predict'
include { SELECT_BEST_ORF                  } from './modules/select_best_orf'
include { METAEUK_PREDICT                  } from './modules/metaeuk'
include { PSAURON_METAEUK                  } from './modules/psauron_metaeuk'
include { GMST_PREDICT                     } from './modules/gmst'
include { PSAURON_GMST                     } from './modules/psauron_gmst'
include { MERGE_PREDICTIONS                } from './modules/merge_predictions'
include { MMSEQS2_CLUSTER_PROTEIN           } from './modules/mmseqs2_cluster_protein'
include { VALIDATE_IDS                     } from './modules/validate_ids'
include { BUSCO as BUSCO_TRINITY           } from './modules/busco'
include { BUSCO as BUSCO_QC                } from './modules/busco'
include { TRANSANNOT                       } from './modules/transannot'
include { THINNING_REPORT                  } from './modules/thinning_report'

// --- Main workflow ---

workflow {

    // Validate required parameters
    Utils.validateParams(params)

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

    // -- Representative selection + taxonomy filter + nucleotide dedup --

    SELECT_REP(ch_trinity, CORSET.out.clust, params.species_label)

    // Taxonomy filter first: fewer sequences for dedup, avoids losing a
    // plant gene because its 95%-identical contaminant was the cluster rep.
    MMSEQS2_TAXONOMY(
        SELECT_REP.out.fasta,
        params.mmseqs2_taxonomy_db,
        params.mmseqs2_search_sens,
        params.filter_taxon,
        params.species_label
    )

    MMSEQS2_CLUSTER(MMSEQS2_TAXONOMY.out.fasta, params.species_label)

    // -- Frameshift correction --

    DIAMOND_BLASTX(MMSEQS2_CLUSTER.out.fasta, params.diamond_db, params.species_label)
    CORRECT_FRAMESHIFTS(MMSEQS2_CLUSTER.out.fasta, DIAMOND_BLASTX.out.tsv, params.species_label)

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

    // -- TD2 best ORF selection --

    SELECT_BEST_ORF(
        TD2_PREDICT.out.psauron_scores,
        TD2_PREDICT.out.pep,
        TD2_PREDICT.out.gff3,
        params.species_label
    )

    // -- MetaEuk ORF prediction (parallel to TD2 branch) --

    METAEUK_PREDICT(
        CORRECT_FRAMESHIFTS.out.fasta,
        params.mmseqs2_swissprot,
        params.species_label
    )

    PSAURON_METAEUK(METAEUK_PREDICT.out.faa, params.species_label)

    // -- GeneMarkS-T ORF prediction (parallel to TD2 + MetaEuk) --

    GMST_PREDICT(CORRECT_FRAMESHIFTS.out.fasta, params.species_label)
    PSAURON_GMST(GMST_PREDICT.out.fnn, params.species_label)

    // -- Merge TD2 + MetaEuk + GeneMarkS-T predictions --

    MERGE_PREDICTIONS(
        SELECT_BEST_ORF.out.faa,
        SELECT_BEST_ORF.out.map,
        METAEUK_PREDICT.out.faa,
        METAEUK_PREDICT.out.map,
        PSAURON_METAEUK.out.scores,
        GMST_PREDICT.out.faa,
        GMST_PREDICT.out.map,
        PSAURON_GMST.out.scores,
        params.min_psauron,
        params.species_label
    )

    // -- Protein-level dedup (95% aa identity) --

    MMSEQS2_CLUSTER_PROTEIN(MERGE_PREDICTIONS.out.faa, params.species_label)

    // -- Gene-level Salmon quant on representative transcripts --

    SALMON_INDEX_FINAL(CORRECT_FRAMESHIFTS.out.fasta)

    SALMON_QUANT_FINAL(
        ch_reads_for_salmon,
        SALMON_INDEX_FINAL.out.index.first(),
        'gene_quant'
    )

    // -- ID validation --

    ch_first_quant_sf = SALMON_QUANT_FINAL.out.quant_dir
        .first()
        .map { qdir -> file("${qdir}/quant.sf") }

    VALIDATE_IDS(
        ch_first_quant_sf,
        MMSEQS2_CLUSTER_PROTEIN.out.faa,
        params.species_label
    )

    // -- BUSCO QC + TransAnnot (parallel, on deduped proteins) --

    BUSCO_QC(MMSEQS2_CLUSTER_PROTEIN.out.faa, params.species_label, 'final')

    TRANSANNOT(
        MMSEQS2_CLUSTER_PROTEIN.out.faa,
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
        MERGE_PREDICTIONS.out.map,
        MMSEQS2_CLUSTER_PROTEIN.out.faa,
        SALMON_QUANT_INITIAL.out.quant_dir.collect(),
        SALMON_QUANT_FINAL.out.quant_dir.collect(),
        BUSCO_TRINITY.out.summary,
        BUSCO_QC.out.summary,
        VALIDATE_IDS.out.report,
        SORTMERNA.out.log.map { _sample_id, logfile -> logfile }.collect(),
        MMSEQS2_TAXONOMY.out.breakdown,
        TRANSANNOT.out.annotation,
        MERGE_PREDICTIONS.out.stats,
        MMSEQS2_CLUSTER.out.stats,
        MMSEQS2_CLUSTER_PROTEIN.out.stats,
        params.species_label
    )

    // -- Published outputs (structure matches README § Output) --

}

// Output publishing is handled by publishDir in conf/base.config.

workflow.onComplete {
    log.info ""
    log.info "Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration          : ${workflow.duration}"
    log.info "Output dir        : ${params.outdir}"
    log.info ""
}
