# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Nextflow DSL2 pipeline (`nf-denovoslim`) that collapses ~1.2M fragmented Trinity de novo transcriptome transcripts into a non-redundant gene set: SuperTranscripts, one-best-protein FASTA, gene-level Salmon quantification, and functional annotation. Runs on the NMBU Orion HPC cluster via SLURM + Apptainer.

Three grass species (BMAX, BMED, FPRA), each with 40 paired-end RNA-seq samples across 10 conditions (5 timepoints x 2 tissues).

## Common Commands

```bash
# Environment setup (required before any Nextflow command)
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate Nextflow
module load Java Anaconda3 singularity

# Submit a species run (14-day SLURM head job)
sbatch run_BMAX.sh
sbatch run_BMED.sh
sbatch run_FPRA.sh

# Dry-run to check syntax without executing
nextflow run main.nf -profile apptainer,slurm -preview

# Check pipeline status / resume a failed run
cd ~/AnnualPerennial/nf-denovoslim/runs/BMAX
nextflow log                           # list previous runs
nextflow run $HOME/AnnualPerennial/nf-denovoslim/main.nf \
    -profile apptainer,slurm -resume   # resume from last checkpoint

# Monitor running SLURM jobs
squeue -u $USER
sacct -j <JOBID> --format=JobID,Elapsed,MaxRSS,State

# Validate Nextflow config (prints resolved params)
nextflow config main.nf -profile apptainer,slurm
```

Work directories are on project storage (not home — quota issues):
`$PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/{BMAX,BMED,FPRA}/work/`

## Architecture

### Pipeline Flow (16 steps)

```
Trinity.fasta ──► MMSEQS2_CLUSTER_NT (97% dedup)
Reads ──► SORTMERNA (rRNA removal)
    │                     │
    └──► SALMON_INDEX/QUANT_INITIAL (--hardFilter --dumpEq)
                          │
                    CORSET (clustering) ──► LACE (SuperTranscripts)
                          │
              MMSEQS2_TAXONOMY_CHUNKED (plant filter)
                          │
              FRAMESHIFT_CORRECTION_CHUNKED
                    ┌─────┴─────┐
              TD2_LONGORFS    SALMON_INDEX/QUANT_FINAL
                    │
          MMSEQS2_SEARCH_CHUNKED (SwissProt + Pfam)
                    │
              TD2_PREDICT ──► SELECT_BEST_ORF
                    │              │
               BUSCO_QC      TRANSANNOT
                    └──────┬───────┘
                   VALIDATE_IDS ──► THINNING_REPORT
```

### Code Organization

- `main.nf` — Workflow logic, all channel wiring, process invocations
- `modules/*.nf` — Individual process definitions (one tool per file)
- `subworkflows/*.nf` — Split→process→merge orchestration for chunked steps
- `conf/base.config` — Per-process CPU/memory/time allocations and retry logic
- `nextflow.config` — Params, container assignments, profiles (apptainer, slurm, docker)
- `bin/` — Python scripts (auto-added to PATH inside containers):
  - `select_best_orf.py` — picks one protein per gene using PSAURON coding-potential scores
  - `correct_frameshifts.py` — applies Diamond blastx-guided indel corrections to SuperTranscripts
  - `thinning_report.py` — generates the final pipeline summary statistics
  - `process_sample_sheet.py` — samplesheet parsing helper
- `containers/` — Custom Dockerfiles and pre-built `.sif` files for TD2, Lace, and Diamond+Python
- `run_{BMAX,BMED,FPRA}.sh` — Per-species SLURM sbatch scripts; each launches from `runs/{SPECIES}/` with work dir on `$PROJECTS`

### DSL2 Module Pattern

Processes are defined generically and aliased in `main.nf` for reuse:
```groovy
include { SALMON_INDEX as SALMON_INDEX_INITIAL } from './modules/salmon_index'
include { SALMON_INDEX as SALMON_INDEX_FINAL   } from './modules/salmon_index'
```
Per-alias configuration (flags, publishDir) is set in `conf/base.config` via `withName:`.

The same pattern applies to `MMSEQS2_SEARCH_CHUNKED`, aliased twice for SwissProt and Pfam searches that run in parallel on the same ORFs:
```groovy
include { MMSEQS2_SEARCH_CHUNKED as MMSEQS2_SEARCH_CHUNKED_SWISSPROT } from './subworkflows/mmseqs2_search_chunked'
include { MMSEQS2_SEARCH_CHUNKED as MMSEQS2_SEARCH_CHUNKED_PFAM      } from './subworkflows/mmseqs2_search_chunked'
```

### Chunked Subworkflow Pattern

Three bottleneck steps (taxonomy, search, frameshift) use identical split→process→merge subworkflows:
1. `modules/X_chunked.nf` defines SPLIT, PROCESS_CHUNK, and MERGE processes
2. `subworkflows/X_chunked.nf` wires them: split FASTA → `.flatten()` → parallel process → `.collect()` → merge
3. Chunk sizes configured in `nextflow.config` params

## Critical DSL2 Patterns (Bug Prevention)

### `.first()` for value channels
When a single-output process feeds a per-sample process, convert with `.first()` or the downstream process runs only once:
```groovy
// CORRECT
SALMON_QUANT_INITIAL(ch_reads, SALMON_INDEX_INITIAL.out.index.first(), 'quant')
// WRONG — runs for only 1 sample due to queue channel min-cardinality
SALMON_QUANT_INITIAL(ch_reads, SALMON_INDEX_INITIAL.out.index, 'quant')
```

### Queue channel fork avoidance
A queue channel consumed by two `.map` operators splits items between them (not duplicated). The pipeline avoids this by:
- Deriving `ch_sample_conditions` for Corset from a fresh `Channel.fromPath(params.samplesheet)` read, not from the same channel used for reads
- Re-deriving condition metadata via `extractCondition()` helper after SortMeRNA, instead of joining back

### Container + DB paths
Large databases (SwissProt, Pfam, TrEMBL) are passed as `val` inputs (path strings), not `path` inputs. This prevents Nextflow from staging multi-GB database files into work directories. The container bind mounts (`-B /mnt/project -B /net/fs-2/scale`) make them accessible.

Note: `eggnog_annotations` is a separate TSV file (`e7_as_e5_annotations.tsv`) passed to TRANSANNOT, distinct from the MMseqs2 eggNOG profile DB (`mmseqs2_eggnog`). Both paths are configured in `nextflow.config`.

## Resource Configuration

All in `conf/base.config`. Memory scales with `task.attempt` for automatic retry on OOM:
```groovy
memory = { check_max( 64.GB * task.attempt, 'memory' ) }
```

Error strategy: retry on OOM/SLURM kill exit codes (104, 134, 137, 139, 140, 143, 247), finish on others. `maxRetries = 2`.

Heaviest processes:
- `MMSEQS2_TAXONOMY_CHUNK`: 32 CPUs, 150-250 GB, 18h
- `MMSEQS2_SEARCH_CHUNK`: 16 CPUs, 80-160 GB, 2h
- `SORTMERNA`: 32 CPUs, 64 GB, maxForks=15

## ID Consistency Chain

Pipeline correctness depends on matching IDs across outputs:
```
SuperTranscript FASTA headers = Salmon quant.sf Name column
                              = species_X.faa headers
                              = TransAnnot queryID column
                              = orf_to_gene_map.tsv gene_id column
```
Enforced by `VALIDATE_IDS` process. Any process that renames or filters sequences must preserve this chain.

## Orion HPC Filesystems

All shared storage is IBM Storage Scale (GPFS) exported via NFS4. Each compute node also has a local XFS disk.

| Volume | Env Var | Path | Backed Up | Purge Policy | Quota (martpali / fjellheimlab) |
|--------|---------|------|-----------|-------------|-------------------------------|
| **Home** | `$HOME` | `/net/fs-2/scale/OrionStore/Home/martpali` | Daily + weekly | None while active | 200G soft / 300G hard per user |
| **Projects** | `$PROJECTS` | `/mnt/project` → OrionStore/Projects (610T) | Daily + weekly | None during project | 40T group quota |
| **Scratch** | `$SCRATCH` | `/mnt/SCRATCH/martpali` → OrionStore/Scratch (50T) | Daily + weekly | **6 months inactive → purged Thursdays** | 500G soft / 1T hard per user |
| **ScratchProjects** | `$SCRATCH_PROJECTS` | `/mnt/ScratchProjects` (200T) | Yes | 6 months | Group quota |
| **Labfiles** | `$LABFILES` | OrionStore/Labfiles | Daily + weekly | None during project | 18T group quota |
| **Work (login)** | `$TMPDIR` | `/work/users/martpali` — local XFS, 187G | **No** | "Every 2 weeks" | None |
| **Work (compute)** | `$TMPDIR` | `/work/users/{user}_{jobid}` — local XFS, **3.5T per node** | **No** | Deleted by SLURM epilog | None |

Backups accessible via `.snapshot` directory under each volume. Check quotas with `myquota -u martpali` or `myquota -g fjellheimlab`.

### $TMPDIR on Compute Nodes

SLURM prolog creates a per-job directory on the node's local 3.5T SSD at `/work/users/${USER}_${SLURM_JOB_ID}` and exports both `$TMPDIR` and `$APPTAINER_TMPDIR`. The epilog runs `rm -rf $TMPDIR` on job completion.

**When to use:** Many parallel jobs reading/writing same files, lots of random I/O, lots of temp files. Copy data in, work locally, copy results out.

**When NOT to use:** Jobs needing files visible across multiple nodes, or when the data exceeds available local space.

**Why scratch is disabled in the pipeline** (`process.scratch = false`): When Nextflow uses node-local scratch, it `mv`s results back to the work dir. `mv` preserves file ownership (`martpali:martpali`), which counts against personal quota instead of `fjellheimlab` group quota, causing quota failures. The `afterScript` that runs `chgrp -R fjellheimlab .` fixes this for the NFS work directory but can't retroactively fix moved files.

**Node-local disk can fill up.** Some nodes have been observed at 97% usage due to stale files from other users. The documented "purge every 2 weeks" policy is not consistently enforced. Clean up your own stale files with:
```bash
# On a compute node (via srun)
rm -rf /work/users/martpali/Rtmp* /work/users/martpali/tmp.* /work/users/martpali/nxf-*
# On login node
rm -rf /work/users/martpali/nxf-*
```

### SLURM Configuration

- Partition: `orion` — 6 compute nodes (cn-31 to cn-35, cn-37), 384 CPUs, ~1.5 TB RAM each, shared
- `afterScript` runs `chgrp -R fjellheimlab . && chmod -R g+rwX .` on every task for quota accounting
- Submit rate limit: 30/min, queue size: 20
- Custom containers (TD2, Lace) are local `.sif` files, not pulled from registries
- `$APPTAINER_TMPDIR` is automatically set by SLURM prolog on compute nodes

## MMseqs2 Taxonomy Internals

The `MMSEQS2_TAXONOMY` process runs `mmseqs taxonomy` against the UniProt/TrEMBL database (~252M sequences, ~132 GB on disk). Understanding its internal behavior is important for interpreting logs and tuning resources.

### Internal target-split mechanism

When the target database index exceeds `--split-memory-limit` (set to 255G), MMseqs2 automatically splits the **target DB** into chunks. For TrEMBL this means 5 splits of ~50M sequences each. Each split builds a ~96 GB prefilter k-mer index in memory, searches **all** query sequences against that split, then discards the index and moves to the next split. This is purely internal to MMseqs2 and is distinct from any pipeline-level query chunking.

Log signature: `Target split mode. Searching through 5 splits.`

### Two-phase taxonomy workflow

`mmseqs taxonomy` runs two search phases internally:

1. **ORF filter** (sensitivity `-s 2`): 6-frame translates all SuperTranscripts into ORFs, runs a fast prefilter+alignment against TrEMBL. Only ORFs with hits pass (e.g., 5.09M ORFs → 391K, ~7.7%). This is a coarse filter, not the final taxonomy assignment. LCA is not applied here.

2. **Full taxonomy search** (sensitivity `-s 7`): Searches only the 391K filtered ORFs with `--lca-search 1` enabled and **LCA mode 3** (weighted). Each target split is searched independently, and LCA merges the results across all splits to assign the final taxonomy.

The log line `LCA mode: 0` that appears in the prefilter phase parameters refers to Phase 1 only. LCA mode 3 is active in Phase 2 where actual taxonomy is assigned.

### Memory implications

Each of the 5 target splits builds a ~96 GB k-mer index. Combined with the query ORFs and alignment workspace, the process typically needs 250-400 GB resident memory. The config starts at 300 GB and retries at 400/500 GB.

## Known Gotchas

1. **Lace caps at 50 transcripts per cluster** — "WARNING: will only take first 50 transcripts" is expected for highly fragmented clusters. Adjustable via `--maxTrans` but 50 is fine.
2. **Taxonomy filter removes no-hit sequences too** — `mmseqs filtertaxdb` drops sequences with no taxonomy hit, not just non-plant. This is intentional (removes contamination).
3. **Corset cluster count scales with samples** — 1 sample → ~50K clusters; 40 samples across 10 conditions → ~244K clusters. This is correct condition-aware behavior.
4. **MMseqs2 taxonomy memory** — TrEMBL DB is ~132 GB on disk; each chunk needs ~55-63 GB resident just for the DB.
5. **Sample naming convention** — `{SPECIES}{IndivID}_{Timepoint}_{Tissue}` (e.g., `BMAX56_T4_L`). The `extractCondition()` helper in `main.nf` takes the last two `_`-delimited parts. Changing sample naming breaks Corset condition grouping.
6. **Frameshift correction uses separate containers** — Diamond blastx runs in the Diamond container, then the Python correction script runs in the BioPython container. These were split because the Diamond container lacks Python3. The `diamond_python` container in `containers/` is a legacy build.
