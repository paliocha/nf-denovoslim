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
                               │
                               ▼
                      SORTMERNA_INDEX
                 (build index from 8 rRNA DBs)
                               │
                               ▼
                         SORTMERNA
                    (per sample, remove rRNA)
                               │
                               ▼
                      filtered reads (non-rRNA)
                               │
      trinity_assembly.fasta   │
                │              │
                ┌──────────────┼───────────────┐
                ▼              ▼               ▼
         MMSEQS2_CLUSTER_NT   SALMON_INDEX  (filtered reads)
           (97% nt dedup)        │               │
                │                ▼               │
                │          SALMON_QUANT_INITIAL ◄┘
                │           (per sample, --hardFilter --dumpEq)
                │                │
                │                ▼
                └──────► CORSET ◄────────────┘
                     (hierarchical clustering on eq classes)
                               │
                               ▼
                           LACE
                     (one SuperTranscript per gene)
                               │
                               ▼
                      MMSEQS2_TAXONOMY
              (taxonomy + filtertaxdb → plant only)
                               │
                               ▼
                   FRAMESHIFT_CORRECTION
              (Diamond blastx → fix assembly indels)
                               │
                ┌──────────────┼──────────────────┐
                ▼              ▼                  ▼
         TD2_LONGORFS    SALMON_INDEX_FINAL  (filtered reads)
                │              │                  │
                ▼              ▼                  │
         MMSEQS2_SEARCH  SALMON_QUANT_FINAL ◄────┘
         (vs SwissProt      (gene-level)
          + Pfam)                │
                │                ▼
                ▼          VALIDATE_IDS
         TD2_PREDICT       (ID consistency)
         (--retain-mmseqs-hits)
                │
                ▼
         SELECT_BEST_ORF
         (one protein per gene,
          mapping file output)
                │
                ├──► species_X.faa
                ├──► best_orfs.gff3
                ├──► orf_to_gene_map.tsv
                │
        ┌───────┴───────────┐
        ▼                   ▼
   BUSCO_QC            TRANSANNOT
                  (.faa vs SwissProt
                   + Pfam + eggNOG7
                   in one step)
                        │
                        ▼
                 THINNING_REPORT
                (pipeline summary)
```

| Step | Process | Tool |
|------|---------|------|
| 0 | rRNA filtering | SortMeRNA 4.3.7 |
| 1 | Nucleotide deduplication (97% id) | MMseqs2 |
| 2–3 | Initial quantification | Salmon 1.10.3 |
| 4 | Hierarchical transcript→gene clustering | Corset 1.09 |
| 5 | Build SuperTranscripts per gene | Lace 1.14.1 |
| 5b | Taxonomy filter (keep Streptophyta only) | MMseqs2 taxonomy + filtertaxdb |
| 5c | Frameshift correction (fix assembly indels) | Diamond blastx + Python |
| 6–9 | ORF prediction with homology support | TD2 + MMseqs2 |
| 10 | Best ORF selection (one protein per gene) | Python/BioPython |
| 11 | Gene-level quantification | Salmon 1.10.3 |
| 12 | ID consistency validation | bash/awk |
| 13 | Protein completeness | BUSCO v6 |
| 14 | Functional annotation | TransAnnot (SwissProt + Pfam + eggNOG7) |
| 15 | Summary report | Python |

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04
- [Apptainer](https://apptainer.org/) (or Singularity/Docker)
- Pre-built MMseqs2 databases: SwissProt, Pfam, eggNOG7, UniProt/TrEMBL (for taxonomy)
- Pre-built Diamond database: UniRef90 (for frameshift correction)
- Pre-built TD2 container (`containers/td2/td2_1.0.8.sif`)
- Pre-built Lace container (`containers/lace/lace_1.14.1_nx2.sif`)

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
| `--outdir` | `./results` | Output directory |

## Pre-building databases

```bash
# SwissProt
mmseqs databases UniProtKB/Swiss-Prot SwissProtDB tmp

# Pfam
mmseqs databases Pfam-A.full PfamDB tmp

# eggNOG7 profiles
mmseqs databases eggNOG eggNOG7DB tmp

# UniProt/TrEMBL (for taxonomy classification — broader coverage than SwissProt)
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

The Lace container uses NetworkX 2.x (compatible with Lace's API) and matplotlib 3.5.x (for seaborn-deep style):

```bash
cd containers/lace
./build.sh
```

This builds `lace_1.14.1_nx2.sif` which fixes the NetworkX 3.x incompatibility in the official biocontainer.

## Output

```
results/
├── sortmerna/                            # Filtered reads (non-rRNA)
├── mmseqs2_nt/                           # Deduplicated assembly
├── clustering/
│   ├── corset-clusters.txt               # Transcript→gene mapping (tx2gene)
│   └── corset-counts.txt                 # Gene-level raw counts
├── supertranscripts/
│   └── supertranscripts.fasta            # One SuperTranscript per gene
├── taxonomy/
│   ├── taxRes_lca.tsv                    # MMseqs2 LCA taxonomy assignments
│   ├── taxRes_report                     # MMseqs2 taxonomy report
│   ├── supertranscripts_filtered.fasta   # Plant-only SuperTranscripts
│   └── taxonomy_filter_stats.txt         # Filter statistics
├── frameshift_correction/
│   ├── supertranscripts_corrected.fasta  # Frameshifts fixed
│   └── frameshift_corrections.log        # Correction statistics
├── td2/                                  # ORF predictions
├── proteins/
│   ├── SPECIES.faa                       # One best protein per gene
│   ├── best_orfs.gff3                    # ORF coordinates
│   └── orf_to_gene_map.tsv              # Gene↔ORF↔PSAURON mapping
├── salmon_final/
│   └── {SAMPLE}_st_quant/quant.sf       # Gene-level quantification
├── busco/                                # Completeness assessment
├── transannot/                           # Functional annotation
└── SPECIES_thinning_report.txt           # Pipeline summary
```

Note: initial Salmon quantification (transcript-level, used internally by Corset) is not published to the results directory.

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

The transcript→gene mapping from Corset is available at `results/clustering/corset-clusters.txt` (two columns: `transcript_id`, `cluster_id`). This maps the original deduplicated Trinity transcripts to their gene clusters and can be used for transcript-level import if needed.

## Author

Martin Paliocha — [NMBU](https://www.nmbu.no/)
