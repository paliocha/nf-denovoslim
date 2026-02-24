# CLAUDE.md

Guidance for Claude Code when working with this repository.

## What This Is

Nextflow DSL2 pipeline (`nf-denovoslim`) that collapses a fragmented Trinity de novo transcriptome into a non-redundant gene set: representative transcripts, merged multi-predictor proteins (TD2 + MetaEuk + GeneMarkS-T), gene-level Salmon quantification, and functional annotation. Runs on SLURM + Apptainer (NMBU Orion HPC). Three grass species: BMAX, BMED, FPRA — 40 paired-end samples each (5 timepoints × 2 tissues).

## Common Commands

```bash
# Environment setup
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate Nextflow
module load Java Anaconda3 singularity

# Submit a species run
sbatch run_BMAX.sh

# Dry-run (requires --params for validation)
nextflow run main.nf -profile apptainer,orion -preview \
  --samplesheet BMAX.samplesheet.csv --trinity_fasta /tmp/x \
  --diamond_db /tmp/x --busco_lineage poales_odb12 \
  --mmseqs2_swissprot /tmp/x --mmseqs2_pfam /tmp/x \
  --mmseqs2_eggnog /tmp/x --mmseqs2_taxonomy_db /tmp/x

# Resume a failed run
cd ~/AnnualPerennial/nf-denovoslim/runs/BMAX
nextflow run $HOME/AnnualPerennial/nf-denovoslim/main.nf -profile apptainer,orion -resume

# Validate resolved config
nextflow config main.nf -profile apptainer,orion
```

Work dirs: `$PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/{BMAX,BMED,FPRA}/work/`

## Architecture

### Pipeline Flow (28 processes)

```
Trinity.fasta + Reads
    │
    ├─► BUSCO_TRINITY (parallel, transcriptome mode)
    │
    ├─► SORTMERNA_INDEX + SORTMERNA (rRNA removal)
    │       │
    │       ▼
    ├─► SALMON_INDEX/QUANT_INITIAL (full Trinity, --dumpEq)
    │       │
    │       ▼
    │   CORSET ──► SELECT_REP (longest transcript per cluster)
    │       │
    │       ▼
    │   MMSEQS2_TAXONOMY (plant filter)
    │       │
    │       ▼
    │   MMSEQS2_CLUSTER (95% nt dedup)
    │       │
    │       ▼
    │   DIAMOND_BLASTX ──► CORRECT_FRAMESHIFTS
    │       │
    │       ├─► TD2_LONGORFS ──► MMSEQS2_SEARCH (SwissProt + Pfam)
    │       │       │
    │       │       ▼
    │       │   TD2_PREDICT ──► SELECT_BEST_ORF
    │       │
    │       ├─► METAEUK_PREDICT ──► PSAURON_METAEUK
    │       │
    │       ├─► GMST_PREDICT ──► PSAURON_GMST
    │       │
    │       │       all three ──► MERGE_PREDICTIONS
    │       │                       │
    │       │                       ▼
    │       │               MMSEQS2_CLUSTER_PROTEIN (95% aa dedup)
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

- `main.nf` — Workflow logic, channel wiring. No param declarations (all in nextflow.config). Calls `Utils.validateParams(params)` at entry.
- `nextflow.config` — All params in single `params {}` block. Container images in `def img = [...]` map (not user-facing params). Profiles, process selectors, manifest.
- `conf/base.config` — Per-process CPU/memory/time via `withName:`. TD2 strategy closures. publishDir rules.
- `conf/orion.config` — Orion HPC site config (SLURM, bind mounts, group quota, scratch)
- `lib/Utils.groovy` — `validateParams()` (fail-fast on missing required params), `extractCondition()` (sample name → condition)
- `modules/*.nf` — One process per file (20 module files, 28 process instances via aliasing)
- `bin/` — Python scripts (auto-added to PATH by Nextflow):
  - `select_best_orf.py` — one protein per gene via PSAURON scores + completeness ranking (TD2 branch)
  - `correct_frameshifts.py` — Diamond blastx-guided indel correction
  - `merge_predictions.py` — 3-way merge of TD2 + MetaEuk + GeneMarkS-T (completeness → PSAURON → length ranking)
  - `select_representative.py` — longest transcript per Corset cluster
  - `metaeuk_select_best.py` — best MetaEuk protein per gene (score ranking)
  - `gmst_select_best.py` — best GeneMarkS-T ORF per gene (completeness → length ranking)
  - `thinning_report.py` — pipeline summary statistics (argparse, named flags)

### DSL2 Patterns

**Module aliasing** — generic processes reused with different config:
```groovy
include { SALMON_INDEX as SALMON_INDEX_INITIAL } from './modules/salmon_index'
include { SALMON_INDEX as SALMON_INDEX_FINAL   } from './modules/salmon_index'
```
Per-alias resources and publishDir set in `conf/base.config` via `withName:`.

**`.first()` for value channels** — single-output process feeding per-sample process:
```groovy
SALMON_QUANT_INITIAL(ch_reads, SALMON_INDEX_INITIAL.out.index.first(), 'quant')
```

**Queue channel fork avoidance** — `ch_sample_conditions` for Corset is derived from a fresh `Channel.fromPath(params.samplesheet)` read, not from the reads channel.

**DB paths as `val` not `path`** — prevents Nextflow from staging multi-GB databases into work dirs. Container bind mounts make them accessible.

**Container map, not params** — container images defined as `def img = [...]` in nextflow.config (not inside `params {}`), so they don't pollute `--help` / CLI namespace. Process selectors reference `img.sortmerna`, `img.td2`, etc.

## Resource Configuration

All in `conf/base.config`. Memory scales with `task.attempt`:
```groovy
memory = { [64.GB * task.attempt, params.max_memory as nextflow.util.MemoryUnit].min() }
```

Error strategy: retry on OOM/SLURM-kill exit codes (104, 134, 135, 137, 139, 140, 143, 247), finish on others. `maxRetries = 2`. Global `process.shell = ['/bin/bash', '-euo', 'pipefail']`.

Heaviest processes: TRANSANNOT (700/1000/1200 GB), MMSEQS2_TAXONOMY (700/1000/1200 GB), DIAMOND_BLASTX (700 GB), CORSET (128 GB/48h).

### Three-way ORF Prediction

Three predictors run in parallel on frameshift-corrected representatives:
1. **TD2** (TransDecoder2) — homology-supported, PSAURON-scored. Uses length-scaling (`-m`/`-M`/`-L`) for short transcripts. Best ORF per gene selected by `select_best_orf.py`.
2. **MetaEuk** — profile-based homology search against SwissProt. Best protein per gene selected by `metaeuk_select_best.py` (score ranking).
3. **GeneMarkS-T** — ab initio self-training (single-threaded, 1 CPU). Best ORF per gene selected by `gmst_select_best.py` (completeness → length).

All three are PSAURON-scored and merged by `merge_predictions.py` with ranking: **completeness → PSAURON → length**. Minimum PSAURON threshold `--min_psauron 0.3`.

### TD2 Strategy Presets

Three presets selectable via `--td2_strategy`:

| Strategy | `-m` | `-M` | `-L` | FDR |
|----------|-----:|-----:|-----:|----:|
| `conservative` | 90 | 70 | 0.7 | 0.05 |
| `standard` | 90 | 50 | 0.5 | 0.10 |
| `aggressive` | 90 | 30 | 0.4 | 0.20 |

When `--td2_strategy` is null (default), individual params (`--td2_min_orf_length`, `--td2_abs_min_orf`, `--td2_length_scale`) are used directly.

## ID Consistency

Pipeline correctness depends on matching IDs:
```
Representative FASTA headers = Salmon quant.sf Name = species.faa headers
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

1. **Taxonomy filter keeps Viridiplantae + no-hit** — `filter_taxon=33090` (Viridiplantae) + taxid 0 (no-hit). Non-plant eukaryotes (fungi, nematodes, oomycetes), bacteria, archaea, viruses are all removed. The `awk '$3 != 1'` intermediate step on `filteredTaxResult.index` is **required** — without it, `filtertaxdb` + `createsubdb` passes everything through (MMseqs2 filtertaxdb zeros data but keeps keys)
2. **Corset cluster count scales with samples** — 40 samples → ~244K clusters (condition-aware)
3. **Sample naming** — `{SPECIES}{ID}_{Timepoint}_{Tissue}`. Non-standard names need explicit `condition` column in samplesheet
4. **Frameshift correction uses two containers** — Diamond (blastx) then BioPython (correction script)
5. **No dedup before Corset** — the pipeline deliberately skips nucleotide deduplication; deduping before Corset destroys multi-mapping signal and produces singleton clusters
6. **DIAMOND flags incompatible with `-F` (frameshift mode)** — two flags break `-F 15`:
   - `-g N` (global ranking): Diamond's ungapped-score ranking cannot be slotted into the 3-frame DP extension pipeline.
   - `--iterate`: runs a sequence of sensitivity steps; the final steps (`default`, `sensitive`) use full matrix extension which `-F` does not support.
   Use `--sensitive` directly (not `--iterate --sensitive`). Performance is controlled by `--top 1`.
7. **TD2 `-M` requires `-L`** — setting `-M` (absolute min) without `-L` (length-scale) has no effect. Both must be set together. The logic: accept ORF if `len >= -m` OR (`len >= -M` AND `len/transcript_len >= -L`).
8. **GeneMarkS-T is single-threaded** — `gmst.pl` has no threading. 1 CPU allocated. Runs ~1-2h total on ~150-200K sequences.
9. **Container map invalidates -resume cache** — changing the container assignment mechanism (e.g. `params.x_container` → `img.x`) changes the process hash even if the resolved image string is identical. This breaks `-resume` for all cached tasks.
10. **Protein dedup is post-merge** — `MMSEQS2_CLUSTER_PROTEIN` runs at 95% aa identity after the 3-way merge, removing near-identical proteins from different predictors that survived the per-gene best-selection.
