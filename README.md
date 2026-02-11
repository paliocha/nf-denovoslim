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

## Pipeline Steps

```
Reads ──► SortMeRNA (rRNA removal)
                │
Trinity.fasta   │
    │           │
    ├───► MMseqs2 97% nt dedup
    │           │
    │     Salmon index + quant (--dumpEq)
    │           │
    └───► Grouper (expression-aware gene clustering)
                │
          SuperTranscripts (one seq per gene)
                │
          TD2 ORF prediction (with SwissProt + Pfam homology)
                │
          Select best ORF per gene (PSAURON scoring)
                │
        ┌───────┼───────┐
        │       │       │
   Salmon    BUSCO   TransAnnot
   (final)    QC    (SwissProt +
                     Pfam + eggNOG)
```

| Step | Process | Tool |
|------|---------|------|
| 0 | rRNA filtering | SortMeRNA 4.3.7 |
| 1 | Nucleotide deduplication (97% id) | MMseqs2 |
| 2–3 | Initial quantification | Salmon 1.10.3 |
| 4 | Expression-aware transcript→gene clustering | Grouper |
| 5 | Merge transcripts per gene | Trinity SuperTranscripts |
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
- Pre-built MMseqs2 databases: SwissProt, Pfam, eggNOG7
- Pre-built TD2 container (`containers/td2/td2_1.0.8.sif`)

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
    --outdir /path/to/results
```

### Samplesheet format

CSV with nf-core/rnaseq-compatible format:

```csv
sample,fastq_1,fastq_2,strandedness
SPECIES01_T1_L,/path/to/R1.fq.gz,/path/to/R2.fq.gz,unstranded
```

The `condition` for Grouper clustering is extracted automatically from the sample name as `{Timepoint}_{Tissue}` (e.g., `T1_L`).

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
```

## Building the TD2 container

The TD2 (TransDecoder 2) container must be built from the included Dockerfile:

```bash
apptainer build containers/td2/td2_1.0.8.sif docker-archive://td2.tar
# or from a machine with Docker:
# docker build -t td2:1.0.8 containers/td2/
# docker save td2:1.0.8 -o td2.tar
# apptainer build td2_1.0.8.sif docker-archive://td2.tar
```

## Output

```
results/
├── sortmerna/              # Filtered reads (non-rRNA)
├── mmseqs2_nt/             # Deduplicated assembly
├── salmon_initial/         # Transcript-level quant (for Grouper)
├── grouper/                # Transcript→gene clustering
├── supertranscripts/       # One sequence per gene
├── td2/                    # ORF predictions
├── proteins/
│   ├── SPECIES.faa         # One best protein per gene
│   ├── best_orfs.gff3      # ORF coordinates
│   └── orf_to_gene_map.tsv # Gene↔ORF↔PSAURON mapping
├── salmon_final/           # Gene-level quantification
├── busco/                  # Completeness assessment
├── transannot/             # Functional annotation
└── report/                 # Pipeline summary
```

## Author

Martin Paliocha — [NMBU](https://www.nmbu.no/)
