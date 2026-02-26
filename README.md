# nf-denovoslim

Nextflow DSL2 pipeline that collapses a fragmented Trinity *de novo* transcriptome assembly into a non-redundant gene set with one protein per gene, using three ORF predictors merged by coding-potential scoring.

## Overview

Given a Trinity assembly and paired-end RNA-seq reads, the pipeline produces:

- **Representative transcript FASTA** — longest transcript per Corset gene cluster
- **Protein FASTA** (`.faa`) — one best protein per gene (3-way merge of TD2 + MetaEuk + GeneMarkS-T, deduplicated at 95% aa identity)
- **Gene-level Salmon quantification** — `quant.sf` with IDs matching the `.faa`
- **Functional annotation** — SwissProt + Pfam + eggNOG via TransAnnot
- **Dual BUSCO assessment** — Trinity baseline (transcriptome mode) + final proteins (protein mode)

### Design rationale

The **full, unfiltered Trinity assembly** goes directly into Salmon (with `--dumpEq`, no `--hardFilter`) so that multi-mapping signal across isoforms is preserved for Corset's hierarchical clustering. Deduplicating transcripts *before* Corset destroys this signal and produces mostly singleton clusters.

After Corset clustering, the **longest transcript** per cluster is selected as the representative (SELECT_REP). This avoids the complexity and runtime of SuperTranscript assembly while retaining the most complete sequence per gene for downstream ORF prediction.

### Taxonomy filtering

Step 4 uses MMseqs2 `taxonomy` against UniRef50 to classify every representative transcript and remove non-plant contaminants (bacteria, fungi, oomycetes, metazoa, viruses, etc.), keeping only Viridiplantae (NCBI taxon 33090) and unclassified sequences (taxon 0 — likely novel or species-specific genes with no database hit).

Key parameter choices:

- **`--lca-mode 3` (weighted LCA)** — each ORF's hits are weighted by alignment score, and the LCA is computed across the top-scoring hits using majority voting (`--majority 0.5`). This avoids the artificially deep assignments of naive LCA while being more robust than top-hit-only classification, especially for genes with patchy taxonomic coverage in UniRef50.
- **`--lca-search false` (default)** — the standard multi-step pipeline (extract ORFs → prefilter → align → collect all hits → then compute LCA). The alternative `--lca-search true` computes LCA on-the-fly during alignment, which is faster but less accurate because it cannot consider all hits before deciding. For a filtering step where false negatives mean losing real genes, accuracy matters more than speed.
- **ORF filter (`--orf-filter 1`, sensitivity 2)** — before the main `-s 7.0` search, a fast low-sensitivity pre-screen discards ORFs with no detectable homology, reducing the query set from millions of short ORFs to only those with plausible hits. This dramatically cuts runtime without affecting the final taxonomy assignments.
- **UniRef50 (not UniRef90)** — 60M clusters vs 300M+ sequences. UniRef50 provides sufficient taxonomic resolution for kingdom-level filtering while fitting in ~500 GB RAM and running 3–4× faster than UniRef90.

### Three-way ORF prediction

Three ORF predictors run in parallel on frameshift-corrected representative transcripts:

1. **TD2** (TransDecoder2) — homology-supported ORF prediction. TD2.LongOrfs extracts candidate ORFs with length-scaling for short transcripts, MMseqs2 searches SwissProt + Pfam for homology evidence, then TD2.Predict scores candidates with PSAURON and rescues borderline ORFs with database hits. Best ORF per gene selected by `select_best_orf.py` (PSAURON + completeness ranking).
2. **MetaEuk** — profile-based homology search against SwissProt. Best protein per gene selected by `metaeuk_select_best.py` (score ranking).
3. **GeneMarkS-T** — ab initio self-training gene finder (single-threaded). Best ORF per gene selected by `gmst_select_best.py` (completeness → length ranking).

All three predictors' outputs are PSAURON-scored independently, then merged by `merge_predictions.py` with ranking: **completeness → protein length → PSAURON score**. This prioritises longer proteins to avoid selecting short internal ORFs when a homology-supported full-length prediction exists. Predictions below `--min_psauron 0.3` are discarded. Genes with zero expression across all samples are dropped; genes below the PSAURON threshold but with mean TPM ≥ 1 are rescued. The merged protein set is then deduplicated at 95% amino acid identity via MMseqs2 clustering to remove near-identical proteins from different predictors.

### TD2 ORF prediction tuning

In *de novo* transcriptome assemblies the median representative transcript is often short (~620 nt for FPRA). A flat minimum ORF length of 90 aa (270 nt) discards legitimate short proteins from short transcripts. TD2 v1.0.8 provides **length-scaling** parameters:

| Parameter | Flag | Description |
|-----------|------|-------------|
| `-m` | `--min-length` | Minimum ORF length (aa) for long transcripts |
| `-M` | `--absolute-min-length` | Absolute minimum ORF length (aa) for short transcripts |
| `-L` | `--length-scale` | Accept a short ORF if it covers ≥ this fraction of the transcript length |

The decision logic is: **accept ORF if** `len ≥ -m` **OR** (`len ≥ -M` **AND** `len / transcript_len ≥ -L`). Setting `-M` without `-L` has no effect — both must be specified together.

Three TD2 tuning strategies are available via `--td2_strategy`:

| Strategy | `-m` | `-M` | `-L` | FDR | Use case |
|----------|-----:|-----:|-----:|----:|----------|
| **Conservative** | 90 | 70 | 0.7 | 0.05 | Reference-quality assemblies with long N50. Minimises false-positive ORFs at the cost of losing short proteins. |
| **Standard** (default) | 90 | 50 | 0.5 | 0.10 | Fragmented *de novo* assemblies (median ~500–800 nt). Recovers short proteins proportional to transcript length while filtering noise from longer sequences. |
| **Aggressive** | 90 | 30 | 0.4 | 0.20 | Highly fragmented assemblies or when maximising sensitivity matters more than precision. Recovers very short ORFs (≥30 aa / 90 nt) at the cost of more false positives passing to TD2.Predict. |

```bash
# Use a named strategy (overrides individual TD2 params)
nextflow run main.nf --td2_strategy aggressive ...

# Or set individual params directly (when --td2_strategy is not set)
nextflow run main.nf --td2_abs_min_orf 60 --td2_length_scale 0.6 ...
```

### Dual BUSCO assessment

The pipeline runs BUSCO twice to measure deduplication effectiveness:

1. **BUSCO Trinity** (transcriptome mode) — baseline on the raw Trinity assembly, typically showing high completeness but extreme duplication (~90% D) due to multiple isoforms per gene.
2. **BUSCO QC** (protein mode) — on the final deduplicated proteins (one per gene after 3-way merge + 95% aa dedup). Duplicated BUSCOs should drop dramatically (to ~10–30% D) since each Corset cluster yields one representative → one best ORF → one protein. Residual duplication reflects genuine paralogs, not assembly artifacts.

Protein-mode BUSCO is the definitive quality metric because (a) it operates directly on the pipeline's deliverable, (b) it avoids MetaEuk gene prediction artifacts that inflate duplication in transcriptome mode, and (c) the one-protein-per-gene design means any remaining duplication is biologically real.

## Pipeline

28 processes across four parallel routes: **Assembly** follows the main representative selection → ORF prediction chain, **Quantification** forks after frameshift correction to re-quantify representatives with Salmon, **Annotation** branches from the merged protein set to TransAnnot + BUSCO QC. **BUSCO Trinity** runs independently on the raw assembly as a baseline.

| Step | Process | Tool |
|------|---------|------|
| 0 | rRNA filtering | SortMeRNA 4.3.7 |
| 1 | Initial quantification (full Trinity, `--dumpEq`) | Salmon 1.10.3 |
| 2 | Transcript-to-gene clustering | Corset 1.10 ([paliocha/Corset](https://github.com/paliocha/Corset)) |
| 3 | Select longest transcript per cluster | Python/BioPython |
| 4 | Taxonomy filter (keep Viridiplantae) | MMseqs2 taxonomy |
| 4b | Nucleotide dedup (95% nt identity) | MMseqs2 cluster |
| 5 | Frameshift correction | Diamond blastx 2.1.22 + BioPython |
| 6a | ORF prediction (homology-supported) | TD2 + MMseqs2 (SwissProt, Pfam) |
| 6b | ORF prediction (profile-based) | MetaEuk (vs SwissProt) |
| 6c | ORF prediction (ab initio) | GeneMarkS-T |
| 7 | PSAURON scoring (each predictor) | TD2 / PSAURON |
| 8 | 3-way merge (completeness → length → PSAURON) | Python/BioPython |
| 8b | Protein dedup (95% aa identity) | MMseqs2 cluster |
| 9 | Gene-level quantification | Salmon 1.10.3 |
| 10 | Protein completeness | BUSCO v6 (protein mode) |
| 11 | Trinity completeness (parallel) | BUSCO v6 (transcriptome mode) |
| 12 | Functional annotation | TransAnnot 4.0.0 |
| 13 | Summary report | Python (argparse) |

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 25.10.0
- [Apptainer](https://apptainer.org/) (or Singularity / Docker)
- Pre-built MMseqs2 databases: SwissProt, Pfam, eggNOG7, UniRef50 taxonomy
- Diamond database: UniRef50 (for frameshift correction)
- Local containers: `containers/td2/td2_1.0.8.sif`, `containers/corset/corset_1.10.sif`, `containers/transannot/transannot_4.0.0_eggnog7.sif`

All other containers (SortMeRNA, Salmon, MMseqs2, Diamond, BioPython, BUSCO, MetaEuk, GeneMarkS-T) are pulled automatically from Biocontainers/Seqera.

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
    --mmseqs2_taxonomy_db /path/to/UniRef50taxdb \
    --diamond_db /path/to/uniref50.dmnd \
    --busco_lineage poales_odb12 \
    --outdir /path/to/results
```

On NMBU Orion, use `-profile apptainer,orion` — database paths are pre-configured in `conf/orion.config`.

### Samplesheet format

```csv
sample,fastq_1,fastq_2,strandedness,condition
SPECIES01_T1_L,/path/to/R1.fq.gz,/path/to/R2.fq.gz,unstranded,T1_L
```

The `condition` column is optional. If omitted, condition is extracted from the sample name as `{Timepoint}_{Tissue}` (last two `_`-separated fields).

### Profiles

| Profile | Description |
|---------|-------------|
| `apptainer` | Apptainer containers |
| `singularity` | Singularity containers |
| `docker` | Docker containers |
| `slurm` | Generic SLURM scheduling |
| `orion` | NMBU Orion HPC (DB paths, SLURM queue, bind mounts, group quota fix) |
| `highmem` | 64 CPUs / 1200 GB RAM cap |
| `standard` | 16 CPUs / 128 GB RAM cap (default) |
| `test` | 4 CPUs / 16 GB RAM cap |

Combine as needed: `-profile apptainer,orion` or `-profile apptainer,slurm,highmem`.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--trinity_fasta` | required | Trinity assembly FASTA |
| `--samplesheet` | required | CSV samplesheet |
| `--species_label` | `species_X` | Output file prefix |
| `--mmseqs2_swissprot` | required | MMseqs2 SwissProt DB |
| `--mmseqs2_pfam` | required | MMseqs2 Pfam DB |
| `--mmseqs2_eggnog` | required | MMseqs2 eggNOG7 profiles DB |
| `--mmseqs2_taxonomy_db` | required | MMseqs2 UniRef50 taxonomy DB |
| `--diamond_db` | required | Diamond UniRef50 DB |
| `--busco_lineage` | required | BUSCO lineage (e.g. `poales_odb12`) |
| `--filter_taxon` | `33090` | NCBI taxon ID to keep (Viridiplantae) |
| `--mmseqs2_search_sens` | `7.0` | MMseqs2 `-s` sensitivity |
| `--min_psauron` | `0.3` | Minimum PSAURON score for merged predictions |
| `--td2_strategy` | `null` | TD2 strategy preset: `conservative`, `standard`, or `aggressive` (overrides individual TD2 params) |
| `--td2_min_orf_length` | `90` | Min ORF length (aa) for long transcripts (`-m`) |
| `--td2_abs_min_orf` | `50` | Absolute min ORF length (aa) for short transcripts (`-M`) |
| `--td2_length_scale` | `0.5` | Accept short ORF if it covers ≥ this fraction of transcript (`-L`) |
| `--td2_strand_specific` | `true` | TD2 strand-specific mode |
| `--sortmerna_db_dir` | `null` | Pre-downloaded rRNA DB dir |
| `--unix_group` | `null` | Unix group for quota accounting |
| `--orion_exclude_nodes` | `null` | SLURM nodes to exclude (e.g. `cn-37`) |
| `--outdir` | `./results` | Output directory |

## Pre-building databases

```bash
# SwissProt
mmseqs databases UniProtKB/Swiss-Prot SwissProtDB tmp

# Pfam
mmseqs databases Pfam-A.full PfamDB tmp

# eggNOG7 profiles
mmseqs databases eggNOG eggNOG7DB tmp

# UniRef50 taxonomy + Diamond databases (downloads prerequisites, then builds)
sbatch scripts/download_uniref50.sh

# Or build individually (if prerequisites already downloaded):
#   sbatch scripts/build_uniref50_diamond.sh
#   sbatch scripts/build_uniref50_mmseqstaxdb.sh
```

## Building containers

**TD2:**
```bash
docker build -t td2:1.0.8 containers/td2/
docker save td2:1.0.8 -o td2.tar
apptainer build containers/td2/td2_1.0.8.sif docker-archive://td2.tar
```

**Corset:**
```bash
cd containers/corset && bash build.sh
```

**TransAnnot:**
```bash
cd containers/transannot && bash build.sh
```

## Output

```
results/
├── clustering/
│   ├── corset-clusters.txt                # tx2gene mapping
│   └── corset-counts.txt                  # Gene-level raw counts
├── representatives/
│   └── representatives.fasta              # Longest transcript per cluster
├── taxonomy/
│   ├── taxRes_lca.tsv
│   ├── representatives_filtered.fasta
│   ├── taxonomy_filter_stats.txt
│   └── taxonomy_breakdown.tsv
├── frameshift_correction/
│   └── frameshift_stats.txt
├── mmseqs2_search/
│   ├── swissprot_alnRes.m8
│   └── pfam_alnRes.m8
├── proteins/
│   ├── SPECIES.faa                        # Final deduplicated proteins
│   ├── merge_map.tsv                      # Gene → predictor source mapping
│   ├── merge_stats.txt                    # Per-predictor contribution counts
│   ├── protein_dedup_stats.txt            # Dedup summary (before/after counts)
│   └── cluster_stats.txt                  # MMseqs2 protein cluster sizes
├── salmon_final/
│   └── {SAMPLE}_quant/quant.sf            # Gene-level quantification
├── qc/
│   ├── busco_trinity/                     # BUSCO transcriptome mode
│   └── busco_final/                       # BUSCO protein mode
├── transannot/
└── SPECIES_thinning_report.txt
```

Initial Salmon quant (used internally by Corset) and SortMeRNA filtered reads are not published.

## Using output with tximport

```r
library(tximport)
library(DESeq2)

files <- list.files("results/salmon_final", pattern = "quant.sf",
                     recursive = TRUE, full.names = TRUE)
names(files) <- gsub("_quant$", "", basename(dirname(files)))

# Representative = gene, so tx2gene is an identity mapping
tx2gene <- data.frame(
  TXNAME = read.delim(files[1])$Name,
  GENEID = read.delim(files[1])$Name
)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
```

The Corset transcript-to-gene map is at `results/clustering/corset-clusters.txt` if transcript-level import is needed.

## Author

Martin Paliocha — [NMBU](https://www.nmbu.no/)
