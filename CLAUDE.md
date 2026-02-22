# CLAUDE.md

Guidance for Claude Code when working with this repository.

## What This Is

Nextflow DSL2 pipeline (`nf-denovoslim`) that collapses a fragmented Trinity de novo transcriptome into a non-redundant gene set: SuperTranscripts, one-best-protein FASTA, gene-level Salmon quantification, and functional annotation. Runs on SLURM + Apptainer (NMBU Orion HPC). Three grass species: BMAX, BMED, FPRA — 40 paired-end samples each (5 timepoints × 2 tissues).

## Common Commands

```bash
# Environment setup
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate Nextflow
module load Java Anaconda3 singularity

# Submit a species run
sbatch run_species.sh BMAX

# Dry-run
nextflow run main.nf -profile apptainer,orion -preview

# Resume a failed run
cd ~/AnnualPerennial/nf-denovoslim/runs/BMAX
nextflow run $HOME/AnnualPerennial/nf-denovoslim/main.nf -profile apptainer,orion -resume

# Validate resolved config
nextflow config main.nf -profile apptainer,orion
```

Work dirs: `$PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/{BMAX,BMED,FPRA}/work/`

## Architecture

### Pipeline Flow

```
Trinity.fasta + Reads
    │
    ├─► BUSCO_TRINITY (parallel, transcriptome mode)
    │
    ├─► SORTMERNA (rRNA removal)
    │       │
    │       ▼
    ├─► SALMON_INDEX/QUANT_INITIAL (full Trinity, --dumpEq)
    │       │
    │       ▼
    │   CORSET ──► LACE (SuperTranscripts)
    │       │
    │       ▼
    │   MMSEQS2_TAXONOMY (plant filter)
    │       │
    │       ▼
    │   DIAMOND_BLASTX ──► CORRECT_FRAMESHIFTS
    │       │
    │       ├─► TD2_LONGORFS ──► MMSEQS2_SEARCH (SwissProt + Pfam)
    │       │       │
    │       │       ▼
    │       │   TD2_PREDICT ──► SELECT_BEST_ORF
    │       │                       │
    │       │               ┌───────┴────────┐
    │       │               ▼                ▼
    │       │           BUSCO_QC         TRANSANNOT
    │       │
    │       ├─► SALMON_INDEX/QUANT_FINAL (gene-level)
    │       │       │
    │       │       ▼
    │       │   VALIDATE_IDS
    │       │
    │       └───────────┬────────────────────┘
    │                   ▼
    └──────────► THINNING_REPORT
```

### Code Organization

- `main.nf` — Workflow logic, channel wiring, input validation
- `modules/*.nf` — One process per file
- `conf/base.config` — Per-process CPU/memory/time via `withName:`
- `conf/orion.config` — Orion HPC site config (DB paths, SLURM, bind mounts, group quota)
- `nextflow.config` — Params, container assignments, profiles
- `bin/` — Python scripts (auto-added to PATH):
  - `select_best_orf.py` — one protein per gene via PSAURON scores
  - `correct_frameshifts.py` — Diamond blastx-guided indel correction
  - `thinning_report.py` — pipeline summary statistics

### DSL2 Patterns

**Module aliasing** — generic processes reused with different config:
```groovy
include { SALMON_INDEX as SALMON_INDEX_INITIAL } from './modules/salmon_index'
include { SALMON_INDEX as SALMON_INDEX_FINAL   } from './modules/salmon_index'
```
Per-alias flags and publishDir set in `conf/base.config` via `withName:`.

**`.first()` for value channels** — single-output process feeding per-sample process:
```groovy
SALMON_QUANT_INITIAL(ch_reads, SALMON_INDEX_INITIAL.out.index.first(), 'quant')
```

**Queue channel fork avoidance** — `ch_sample_conditions` for Corset is derived from a fresh `Channel.fromPath(params.samplesheet)` read, not from the reads channel.

**DB paths as `val` not `path`** — prevents Nextflow from staging multi-GB databases into work dirs. Container bind mounts make them accessible.

## Resource Configuration

All in `conf/base.config`. Memory scales with `task.attempt`:
```groovy
memory = { check_max( 64.GB * task.attempt, 'memory' ) }
```

Error strategy: retry on OOM/SLURM-kill exit codes (104, 134, 137, 139, 140, 143, 247), finish on others. `maxRetries = 2`. Global `process.shell = ['/bin/bash', '-euo', 'pipefail']`.

Heaviest processes: TRANSANNOT (500/700/1000 GB), MMSEQS2_TAXONOMY (580/700/900 GB), CORSET (128 GB/48h).

## ID Consistency

Pipeline correctness depends on matching IDs:
```
SuperTranscript FASTA headers = Salmon quant.sf Name = species.faa headers
```
Enforced by `VALIDATE_IDS`. Any process that renames/filters sequences must preserve this chain.

## Orion HPC

### Filesystems

| Volume | Path | Notes |
|--------|------|-------|
| Home | `$HOME` | 200G quota, backed up |
| Projects | `$PROJECTS` = `/mnt/project` | Group quota, backed up |
| Scratch | `$SCRATCH` | 500G quota, purged after 6mo inactive |
| Local SSD | `$TMPDIR` | Per-job on compute, `/work/users/$USER` on login |

`$TMPDIR` on compute nodes is created by SLURM per-job and auto-cleaned. Node-local, not visible across nodes. 3.5–11 TB per node.

Pipeline uses `process.scratch = '$TMPDIR'`: Nextflow stages inputs to node-local SSD, runs task, rsync's outputs back to NFS.

### SLURM

- Partition: `orion` — 6 nodes (cn-31 to cn-35, cn-37), 384 CPUs, ~1.5 TB RAM each
- Submit rate: 30/min, queue size: 20
- Exclude nodes: `--orion_exclude_nodes cn-37`

### sbatch wrapper & BASH_ENV

Orion's site sbatch wrapper (`/cluster/software/slurm/site/bin/sbatch`) injects `BASH_ENV` pointing to a Slurm env-sourcing script. The script (`slurm_bash_env.sh`) has three recursion guards built in (`ORION_SLURM_ENV_APPLIED`, `_SLURM_BASH_ENV_SOURCED`, and `unset BASH_ENV` at the end), so no manual `unset BASH_ENV` is needed in run scripts.

## Known Gotchas

1. **Lace caps at 50 transcripts per cluster** — expected for fragmented clusters
2. **Taxonomy filter keeps Viridiplantae + no-hit** — `filter_taxon=33090` (Viridiplantae) + taxid 0 (no-hit). Non-plant eukaryotes (fungi, nematodes, oomycetes), bacteria, archaea, viruses are all removed. The `awk '$3 != 1'` intermediate step on `filteredTaxResult.index` is **required** — without it, `filtertaxdb` + `createsubdb` passes everything through (MMseqs2 filtertaxdb zeros data but keeps keys)
3. **Corset cluster count scales with samples** — 40 samples → ~244K clusters (condition-aware)
4. **Sample naming** — `{SPECIES}{ID}_{Timepoint}_{Tissue}`. Non-standard names need explicit `condition` column in samplesheet
5. **Frameshift correction uses two containers** — Diamond (blastx) then BioPython (correction script)
6. **No dedup before Corset** — the pipeline deliberately skips nucleotide deduplication; deduping before Corset destroys multi-mapping signal and produces singleton clusters
7. **DIAMOND flags incompatible with `-F` (frameshift mode)** — two flags break `-F 15`:
   - `-g N` (global ranking): Diamond's ungapped-score ranking cannot be slotted into the 3-frame DP extension pipeline.
   - `--iterate`: runs a sequence of sensitivity steps; the final steps (`default`, `sensitive`) use full matrix extension which `-F` does not support. Error: `Frameshift alignment does not support full matrix extension`.
   Use `--sensitive` directly (not `--iterate --sensitive`). Performance is controlled by `--top 1` (limits gapped extensions to near-best targets).
