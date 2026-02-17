# nf-denovoslim

A Nextflow DSL2 pipeline that collapses a fragmented Trinity *de novo* transcriptome assembly into a non-redundant gene set.

## Overview

Given a Trinity assembly and paired-end RNA-seq reads, the pipeline produces:

- **SuperTranscript FASTA** — one consensus sequence per gene
- **Protein FASTA** (`.faa`) — one best protein per gene (for OrthoFinder / functional annotation)
- **GFF3** — ORF coordinates on SuperTranscripts
- **Gene-level Salmon quantification** — `quant.sf` with IDs matching the `.faa`
- **Functional annotation** — SwissProt + Pfam + eggNOG via TransAnnot
- **BUSCO completeness** assessment

## Pipeline DAG

```
                        (input reads)
                               |
                               v
                      SORTMERNA_INDEX
                 (build index from 8 rRNA DBs)
                               |
                               v
                         SORTMERNA
                    (per sample, remove rRNA)
                               |
                               v
                      filtered reads (non-rRNA)------------------+
                               |                                 |
      trinity_assembly.fasta   |                                 |
                |              |                                 |
                +--------------+---------------+                 |
                v              v               v                 |
         MMSEQS2_CLUSTER_NT   SALMON_INDEX  (for initial quant)  |
           (97% nt dedup)        |               |               |
                |                v               |               |
                |          SALMON_QUANT_INITIAL <-+              |
                |           (per sample, --hardFilter --dumpEq)  |
                |                |                               |
                |                v                               |
                +--------> CORSET <-----------+                  |
                     (hierarchical clustering on eq classes)     |
                               |                                 |
                               v                                 |
                              LACE                               |
                     (one SuperTranscript per gene)              |
                               |                                 |
                               v                                 |
                  MMSEQS2_TAXONOMY_CHUNKED                       |
          (split -> taxonomy + filter -> merge)                  |
                  (keep Streptophyta only)                       |
                               |                                 |
                               v                                 |
               FRAMESHIFT_CORRECTION_CHUNKED                     |
          (split -> Diamond blastx -> correct -> merge)          |
               (fix assembly indels)                             |
                               |                                 |
                +--------------+------------------+              |
                v              v                  v              |
         TD2_LONGORFS    SALMON_INDEX_FINAL  (for final quant)   |
                |              |                  |              |
                v              v                  v              |
  MMSEQS2_SEARCH_CHUNKED  SALMON_QUANT_FINAL <-------------------+
   (split -> search         (gene-level)
    -> merge)                    |
   (vs SwissProt + Pfam)        v
                |          VALIDATE_IDS
                v          (ID consistency)
         TD2_PREDICT
         (--retain-mmseqs-hits)
                |
                v
         SELECT_BEST_ORF
         (one protein per gene,
          mapping file output)
                |
                +---> species_X.faa
                +---> best_orfs.gff3
                +---> orf_to_gene_map.tsv
                |
        +-------+-----------+
        v                   v
   BUSCO_QC            TRANSANNOT
                  (.faa vs SwissProt
                   + Pfam + eggNOG7
                   in one step)
                        |
                        v
                 THINNING_REPORT
                (pipeline summary)
```

| Step | Process | Tool |
|------|---------|------|
| 0 | rRNA filtering | SortMeRNA 4.3.7 |
| 1 | Nucleotide deduplication (97% id) | MMseqs2 |
| 2-3 | Initial quantification | Salmon 1.10.3 |
| 4 | Hierarchical transcript-to-gene clustering | Corset 1.09 |
| 5 | Build SuperTranscripts per gene | Lace 1.14.1 |
| 5b | Taxonomy filter (keep Streptophyta only, chunked) | MMseqs2 taxonomy + filtertaxdb |
| 5c | Frameshift correction (fix assembly indels, chunked) | Diamond blastx + Python |
| 6-9 | ORF prediction with homology support (chunked search) | TD2 + MMseqs2 |
| 10 | Best ORF selection (one protein per gene) | Python/BioPython |
| 11 | Gene-level quantification | Salmon 1.10.3 |
| 12 | ID consistency validation | bash/awk |
| 13 | Protein completeness | BUSCO v6 |
| 14 | Functional annotation | TransAnnot (SwissProt + Pfam + eggNOG7) |
| 15 | Summary report | Python |

Steps 5b, 5c, and 6-9 use data-level parallelization: input sequences are split into chunks, processed in parallel, and merged. This provides 3-10x wall-clock speedup on HPC clusters.

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04
- [Apptainer](https://apptainer.org/) (or Singularity/Docker)
- Pre-built MMseqs2 databases: SwissProt, Pfam, eggNOG7, UniProt/TrEMBL (for taxonomy)
- Pre-built Diamond database: UniRef90 (for frameshift correction)
- Pre-built TD2 container (`containers/td2/td2_1.0.8.sif`)
- Pre-built Lace container (`containers/lace/lace_1.14.1_patched.sif`)

## Usage

```bash
nextflow run main.nf \
    -profile apptainer,slurm \
    --trinity_fasta /path/to/Trinity.fasta \
    --samplesheet samples.csv \
    --species_label SPECIES \
    --mmseqs2_swissprot /path/to/SwissProtDB \
    --mmseqs2_pfam /path/to/PfamDB \
    --mmseqs2_eggnog /path/to/eggNOG7DB \
    --mmseqs2_taxonomy_db /path/to/UniProtTrEMBLtaxdb \
    --diamond_db /path/to/uniref90.dmnd \
    --outdir /path/to/results
```

### Samplesheet format

CSV with nf-core/rnaseq-compatible format:

```csv
sample,fastq_1,fastq_2,strandedness
SPECIES01_T1_L,/path/to/R1.fq.gz,/path/to/R2.fq.gz,unstranded
```

The `condition` for Corset clustering is extracted automatically from the sample name as `{Timepoint}_{Tissue}` (e.g., `T1_L`).

### Profiles

| Profile | Description |
|---------|-------------|
| `apptainer` | Run with Apptainer containers |
| `singularity` | Run with Singularity containers |
| `slurm` | Submit jobs to SLURM scheduler |
| `docker` | Run with Docker containers |

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--trinity_fasta` | required | Path to Trinity assembly FASTA |
| `--samplesheet` | required | CSV samplesheet (see format above) |
| `--species_label` | `species_X` | Prefix for output files |
| `--mmseqs2_nt_id` | `0.97` | Nucleotide dedup identity threshold |
| `--mmseqs2_nt_cov` | `0.8` | Nucleotide dedup coverage threshold |
| `--filter_taxon` | `35493` | NCBI taxon ID to keep (includes all descendants) |
| `--diamond_db` | required | Path to Diamond database (e.g., UniRef90) for frameshift correction |
| `--busco_lineage` | `poales_odb12` | BUSCO lineage dataset |
| `--sortmerna_db_dir` | `null` | Pre-downloaded rRNA database dir (skips download) |
| `--outdir` | `./results` | Output directory |

#### Chunking parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--search_orf_chunk_size` | `40000` | ORFs per chunk for MMseqs2 search |
| `--taxonomy_chunk_size` | `2000` | SuperTranscripts per chunk for MMseqs2 taxonomy |
| `--max_parallel_search_chunks` | `8` | Max simultaneous search chunk jobs |
| `--max_parallel_taxonomy_chunks` | `8` | Max simultaneous taxonomy chunk jobs |

## Pre-building databases

```bash
# SwissProt
mmseqs databases UniProtKB/Swiss-Prot SwissProtDB tmp

# Pfam
mmseqs databases Pfam-A.full PfamDB tmp

# eggNOG7 profiles
mmseqs databases eggNOG eggNOG7DB tmp

# UniProt/TrEMBL (for taxonomy classification -- broader coverage than SwissProt)
mmseqs databases UniProtKB/TrEMBL UniProtTrEMBLtaxdb tmp

# UniRef90 Diamond database (for frameshift correction)
# Download and build using the provided script:
sbatch build_uniref90_diamond.sh
```

## Building custom containers

### TD2 container

The TD2 (TransDecoder 2) container must be built from the included Dockerfile:

```bash
apptainer build containers/td2/td2_1.0.8.sif docker-archive://td2.tar
# or from a machine with Docker:
# docker build -t td2:1.0.8 containers/td2/
# docker save td2:1.0.8 -o td2.tar
# apptainer build td2_1.0.8.sif docker-archive://td2.tar
```

### Lace container

The Lace container patches Lace 1.14.1 for NetworkX 3.x compatibility (`.node` -> `.nodes`) and pins matplotlib <3.6 (for the `seaborn-deep` style):

```bash
cd containers/lace
./build.sh
```

This builds `lace_1.14.1_patched.sif`.

## Output

```
results/
├── clustering/
│   ├── corset-clusters.txt               # Transcript-to-gene mapping (tx2gene)
│   └── corset-counts.txt                 # Gene-level raw counts
├── supertranscripts/
│   └── supertranscripts.fasta            # One SuperTranscript per gene
├── taxonomy/
│   ├── taxRes_lca.tsv                    # MMseqs2 LCA taxonomy assignments
│   ├── supertranscripts_filtered.fasta   # Plant-only SuperTranscripts
│   └── taxonomy_filter_stats.txt         # Filter statistics
├── frameshift_correction/
│   └── frameshift_stats.txt              # Correction statistics
├── mmseqs2_search/
│   ├── swissprot_alnRes.m8               # SwissProt homology hits
│   └── pfam_alnRes.m8                    # Pfam homology hits
├── proteins/
│   └── SPECIES.faa                       # One best protein per gene
├── annotation/
│   ├── best_orfs.gff3                    # ORF coordinates
│   └── orf_to_gene_map.tsv              # Gene-to-ORF-to-PSAURON mapping
├── salmon_final/
│   └── {SAMPLE}_st_quant/quant.sf       # Gene-level quantification
├── qc/
│   └── busco/                            # BUSCO completeness assessment
├── transannot/                           # Functional annotation
└── SPECIES_thinning_report.txt           # Pipeline summary
```

Note: initial Salmon quantification (transcript-level, used internally by Corset) and SortMeRNA filtered reads are not published to the results directory.

## Using output with tximport (R/DESeq2)

The gene-level Salmon quantification on SuperTranscripts can be imported directly into R for differential expression analysis. Each SuperTranscript corresponds to one gene, so the `quant.sf` Name column matches the protein IDs in the `.faa`.

```r
library(tximport)
library(DESeq2)

# List all gene-level quant.sf files
files <- list.files("results/salmon_final", pattern = "quant.sf",
                     recursive = TRUE, full.names = TRUE)
names(files) <- gsub("_st_quant$", "", basename(dirname(files)))

# tx2gene: identity mapping (SuperTranscript = gene)
tx2gene <- data.frame(
  TXNAME = read.delim(files[1])$Name,
  GENEID = read.delim(files[1])$Name
)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
```

The transcript-to-gene mapping from Corset is available at `results/clustering/corset-clusters.txt` (two columns: `transcript_id`, `cluster_id`). This maps the original deduplicated Trinity transcripts to their gene clusters and can be used for transcript-level import if needed.

## Author

Martin Paliocha -- [NMBU](https://www.nmbu.no/)
