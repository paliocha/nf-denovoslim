# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Nextflow DSL2 pipeline (`nf-denovoslim`) that collapses ~1.2M fragmented Trinity de novo transcriptome transcripts into a non-redundant gene set: SuperTranscripts, one-best-protein FASTA, gene-level Salmon quantification, and functional annotation. Designed for HPC clusters with SLURM + Apptainer, with site-specific config isolated in `conf/orion.config`.

Three grass species (BMAX, BMED, FPRA), each with 40 paired-end RNA-seq samples across 10 conditions (5 timepoints x 2 tissues).

## Common Commands

```bash
# Environment setup (required before any Nextflow command)
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate Nextflow
module load Java Anaconda3 singularity

# Submit a species run (14-day SLURM head job)
sbatch run_species.sh BMAX
sbatch run_species.sh BMED
sbatch run_species.sh FPRA

# Dry-run to check syntax without executing
nextflow run main.nf -profile apptainer,orion -preview

# Check pipeline status / resume a failed run
cd ~/AnnualPerennial/nf-denovoslim/runs/BMAX
nextflow log                           # list previous runs
nextflow run $HOME/AnnualPerennial/nf-denovoslim/main.nf \
    -profile apptainer,orion -resume   # resume from last checkpoint

# Monitor running SLURM jobs
squeue -u $USER
sacct -j <JOBID> --format=JobID,Elapsed,MaxRSS,State

# Validate Nextflow config (prints resolved params)
nextflow config main.nf -profile apptainer,orion
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

- `main.nf` — Workflow logic, all channel wiring, process invocations, input validation
- `modules/*.nf` — Individual process definitions (one tool per file)
- `subworkflows/*.nf` — Split→process→merge orchestration for chunked steps
- `conf/base.config` — Per-process CPU/memory/time allocations and retry logic
- `conf/orion.config` — NMBU Orion HPC site-specific settings (DB paths, SLURM queue, filesystem mounts, group quota fix)
- `nextflow.config` — Params (all DB defaults `null`), container assignments, profiles (apptainer, slurm, orion, docker, highmem, standard, test)
- `bin/` — Python scripts (auto-added to PATH inside containers):
  - `select_best_orf.py` — picks one protein per gene using PSAURON coding-potential scores
  - `correct_frameshifts.py` — applies Diamond blastx-guided indel corrections to SuperTranscripts
  - `thinning_report.py` — generates the final pipeline summary statistics
  - `process_sample_sheet.py` — samplesheet parsing helper
- `containers/` — Custom Dockerfiles and pre-built `.sif` files for TD2, Lace, and Diamond+Python
- `run_species.sh` — Consolidated SLURM sbatch script; takes species code as CLI arg (e.g., `sbatch run_species.sh BMAX`)
- `run_{BMAX,BMED,FPRA}.sh` — **Deprecated** per-species SLURM scripts (use `run_species.sh` instead)

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
Large databases (SwissProt, Pfam, UniRef90) are passed as `val` inputs (path strings), not `path` inputs. This prevents Nextflow from staging multi-GB database files into work directories. The Orion container bind mounts (`-B /mnt/project -B /net/fs-2/scale`) in `conf/orion.config` make them accessible.

All DB path params default to `null` in `nextflow.config` and are set either by the `orion` profile (`conf/orion.config`) or via CLI flags. The pipeline validates that all required DB params are non-null at startup.

Note: `eggnog_annotations` is a separate TSV file (`e7_as_e5_annotations.tsv`) passed to TRANSANNOT, distinct from the MMseqs2 eggNOG profile DB (`mmseqs2_eggnog`). Both paths must be provided.

## Resource Configuration

All in `conf/base.config`. Memory scales with `task.attempt` for automatic retry on OOM:
```groovy
memory = { check_max( 64.GB * task.attempt, 'memory' ) }
```

Error strategy: retry on OOM/SLURM kill exit codes (104, 134, 137, 139, 140, 143, 247), finish on others. `maxRetries = 2`.

Resource caps are set by profiles: `standard` (16 CPUs, 128 GB — default), `highmem` (32 CPUs, 1200 GB), or `test` (4 CPUs, 16 GB). The `orion` profile includes its own caps (32 CPUs, 1200 GB, 14 days).

Global bash strict mode is enabled: `process.shell = ['/bin/bash', '-euo', 'pipefail']`.

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

All shared storage is IBM Storage Scale (GPFS) exported via NFS4. Each compute node also has a local XFS disk mounted at `/work` (3.5–11 TB depending on node).

| Volume | Env Var | Path | Backed Up | Purge Policy | Quota |
|--------|---------|------|-----------|-------------|-------|
| **Home** | `$HOME` | `/net/fs-2/scale/OrionStore/Home/$USER` | Yes | None while active | 200G soft / 300G hard per user |
| **Projects** | `$PROJECTS` | `/mnt/project` → OrionStore/Projects | Yes | None during project | Group quota |
| **Scratch** | `$SCRATCH` | `/mnt/SCRATCH/$USER` | Yes | **6 months inactive → purged** | 500G soft / 1T hard per user |
| **ScratchProjects** | `$SCRATCH_PROJECTS` | `/mnt/ScratchProjects` | Yes | 6 months | Group quota |
| **Labfiles** | `$LABFILES` | OrionStore/Labfiles | Yes | None during project | Group quota |
| **Local work** | `$TMPDIR` | `/work/users/$USER` — node-local XFS | **No** | Not automatically cleaned | None |

### $TMPDIR — node-local SSD

`$TMPDIR` resolves to `/work/users/$USER` on **both** login and compute nodes. The path is identical but each node has its own physical disk — files on the login node's `/work` are NOT visible from compute nodes and vice versa.

**Key facts:**
- **Not per-job** — SLURM does NOT create or clean up per-job subdirectories. `$TMPDIR` is a persistent user directory.
- **No automatic cleanup** — stale files accumulate unless you clean them yourself. The documented "purge every 2 weeks" on login nodes is not consistently enforced.
- **3.5–11 TB per node** — capacity varies; some nodes have been observed at 97% usage due to other users' stale files.
- **When to use:** Many parallel jobs reading/writing the same data, lots of random I/O, lots of temporary files.
- **When NOT to use:** Jobs needing files visible across multiple nodes, or data exceeding available space.

**How scratch works in this pipeline** (`process.scratch = '$TMPDIR'`): Nextflow creates a unique temp directory under `$TMPDIR` (e.g., `/work/users/martpali/nxf-XXXXXX`) for each task. It stages inputs there, runs the task, then `rsync`s declared outputs back to the NFS work directory. On success, Nextflow removes the temp dir. On failure, the temp dir persists for debugging.

The `beforeScript` sets setgid+umask on the scratch dir so task-created files inherit the project group, and also sets setgid on the NFS work dir via `$NXF_TASK_WORKDIR` so files rsync'd back inherit the correct group. Processes that copy large DB files (MMSEQS2_TAXONOMY, DIAMOND_BLASTX) do so to CWD explicitly since DBs are passed as `val` path strings, not `path` inputs.

**Cleanup of stale Nextflow scratch dirs:**
```bash
# On a compute node (via srun)
srun --partition=orion --time=5:00 --ntasks=1 bash -c 'rm -rf /work/users/$USER/nxf-*'
# On login node
rm -rf /work/users/$USER/nxf-*
```

### SLURM Configuration

- All Orion-specific settings are in `conf/orion.config` (loaded by `-profile orion`).
- Partition: `orion` — 6 compute nodes (cn-31 to cn-35, cn-37), 384 CPUs, ~1.5 TB RAM each, shared
- `scratch = '$TMPDIR'` — all tasks run on node-local SSD (`/work/users/$USER`); Nextflow creates per-task `nxf-XXXXXX` dirs and rsync's outputs back to NFS work dir
- `beforeScript` sets setgid+umask on scratch dir AND sets setgid on NFS work dir (via `$NXF_TASK_WORKDIR`) for correct group ownership; `afterScript` runs `chgrp -R <group> .` as safety net — both conditional on `params.unix_group`
- Submit rate limit: 30/min, queue size: 20
- Problem nodes can be excluded via `--orion_exclude_nodes cn-37`
- Custom containers (TD2, Lace) are local `.sif` files, not pulled from registries
- Container bind mounts include `-B /work:/work` so scratch dirs are accessible inside Apptainer

## MMseqs2 Taxonomy Internals

The `MMSEQS2_TAXONOMY` process runs `mmseqs taxonomy` against the UniRef90 database (~90M representative sequences, ~25 GB on disk). UniRef90 replaced the full TrEMBL database (~252M sequences) to eliminate target-splitting and dramatically reduce runtime.

### Why UniRef90 instead of full TrEMBL

The taxonomy step is used only as a coarse kingdom-level filter (Viridiplantae, taxon 35493). At 90% identity clustering, taxonomic lineages are virtually identical within each cluster, so there is no loss of filter accuracy. UniRef90 gives a ~3.5× smaller k-mer index that fits in RAM without any target splitting.

### Memory formula (from MMseqs2 wiki)

```
M = (7 * N * L + 8 * a^k) bytes
```

For UniRef90 (~90M seqs, avg length ~330): M ≈ 208 GB + 10 GB = ~218 GB.
With `--split-memory-limit` at 85% of task memory, 400 GB allocation → 340 GB limit → no splits.

### Two-phase taxonomy workflow

`mmseqs taxonomy` runs two search phases internally:

1. **ORF filter** (sensitivity `-s 2`): 6-frame translates all SuperTranscripts into ORFs, runs a fast prefilter+alignment against the DB. Only ORFs with hits pass (~7-8%). This is a coarse filter, not the final taxonomy assignment.

2. **Full taxonomy search** (sensitivity `-s 7`): Searches only the filtered ORFs with `--lca-search 1` enabled and **LCA mode 3** (weighted). LCA merges results to assign the final taxonomy.

### Previous setup (TrEMBL, for reference)

The original TrEMBL DB had 252M sequences. Its k-mer index was ~554 GB, requiring 5 target splits at 255G `--split-memory-limit`. This introduced a massive result-merging I/O overhead and 16+ hour runtimes. The switch to UniRef90 was made in Feb 2026.

## Known Gotchas

1. **Lace caps at 50 transcripts per cluster** — "WARNING: will only take first 50 transcripts" is expected for highly fragmented clusters. Adjustable via `--maxTrans` but 50 is fine.
2. **Taxonomy filter removes no-hit sequences too** — `mmseqs filtertaxdb` drops sequences with no taxonomy hit, not just non-plant. This is intentional (removes contamination).
3. **Corset cluster count scales with samples** — 1 sample → ~50K clusters; 40 samples across 10 conditions → ~244K clusters. This is correct condition-aware behavior.
4. **MMseqs2 taxonomy memory** — UniRef90 DB k-mer index is ~218 GB; needs ~400 GB allocation to avoid target splitting. Config retry ladder: 150/200/250 GB for chunked mode.
5. **Sample naming convention** — `{SPECIES}{IndivID}_{Timepoint}_{Tissue}` (e.g., `BMAX56_T4_L`). The `extractCondition()` helper in `main.nf` takes the last two `_`-delimited parts. If sample names don't follow this convention, add a `condition` column to the samplesheet CSV to provide explicit Corset groupings.
6. **Frameshift correction uses separate containers** — Diamond blastx runs in the Diamond container, then the Python correction script runs in the BioPython container. These were split because the Diamond container lacks Python3. The `diamond_python` container in `containers/` is a legacy build.
