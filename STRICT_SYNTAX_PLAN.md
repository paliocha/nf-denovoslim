# Strict Syntax / Typed Params / Typed Processes — Refactoring Plan

## Executive Summary

Nextflow 25.10 introduced three progressive language improvements:

| Feature | Maturity | Requires v2 parser? | Our readiness |
|---------|----------|---------------------|---------------|
| **Strict syntax** | Stable (default in 26.04) | Yes (`NXF_SYNTAX_PARSER=v2`) | 95 % — only 2 issues |
| **Typed params** | Stable (25.10+) | Yes | Ready after strict syntax |
| **Typed processes** | **Preview** (`nextflow.preview.types`) | Yes | Moderate rewrite |

Adopting strict syntax is low-risk and already nearly done. Typed params is a
natural follow-up. Typed processes is a larger commitment that should wait until
the feature exits preview.

---

## Current State — Audit Results

The codebase was audited against every strict-syntax rule. **Two issues** remain;
everything else passes cleanly.

### Issues to Fix

| # | File | Lines | Issue | Severity |
|---|------|-------|-------|----------|
| 1 | `main.nf` | 280–286 | `workflow.onComplete` is a **top-level statement** outside the `workflow {}` block | Must-fix for v2 |
| 2 | `main.nf` | 38–44 | Top-level `def extractCondition()` — mixing declarations with statements | Should-fix (deprecated) |

### Already Compliant

- All 17 processes have explicit `script:` labels ✅
- No implicit `it` closures (all use named params: `row ->`, `sc ->`, etc.) ✅
- No `for`/`while` loops, `switch`, `import`, class declarations, spread operators ✅
- No `shell:` process sections — all use `script:` ✅
- All lowercase `channel` (never `Channel`) ✅
- No `addParams` in any include ✅
- No `when:` sections in processes ✅
- No `params.*` references inside process bodies ✅
- All config files use standard patterns ✅

---

## Phase 1 — Strict Syntax (low-risk, do now)

**Goal:** Enable `NXF_SYNTAX_PARSER=v2` and pass `nextflow lint`.

### Step 1.1: Move `workflow.onComplete` inside the workflow block

The v2 parser forbids mixing top-level statements with declarations. Move the
handler inside `workflow {}` using the new `onComplete:` section syntax (25.10+):

```diff
 workflow {
     main:
     // ...existing pipeline logic...

+    onComplete:
+    log.info ""
+    log.info "Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
+    log.info "Duration          : ${workflow.duration}"
+    log.info "Output dir        : ${params.outdir}"
+    log.info ""

     publish:
     // ...existing publish assignments...
 }

-workflow.onComplete {
-    log.info ""
-    log.info "Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
-    ...
-}
```

> **Note:** The `onComplete:` section is a new 25.10 feature that replaces
> `workflow.onComplete = { ... }`.

### Step 1.2: Move `extractCondition()` to `lib/` or inline

Option A (**recommended**): Move to a Groovy helper class in `lib/`:

```
lib/
└── Utils.groovy
```

```groovy
class Utils {
    static String extractCondition(String sample_name) {
        def parts = sample_name.split('_')
        if (parts.size() < 3) return sample_name
        return "${parts[-2]}_${parts[-1]}"
    }
}
```

Then in `main.nf`, replace `extractCondition(...)` with `Utils.extractCondition(...)`.

Option B: Keep it as a top-level `def` in `main.nf` — this still works under
v2 because the declaration is a function (not a statement), but the language
server will flag it as "should be in lib/".

### Step 1.3: Enable the v2 parser

In run scripts (`run_BMAX.sh`, etc.) or in `nextflow.config`:

```bash
export NXF_SYNTAX_PARSER=v2
```

Or add to `nextflow.config`:

```groovy
env {
    NXF_SYNTAX_PARSER = 'v2'
}
```

### Step 1.4: Validate

```bash
nextflow lint main.nf
nextflow run main.nf -preview   # dry-run
```

### Estimated effort: ~30 minutes, zero risk to running pipelines.

---

## Phase 2 — Typed Params (low-risk, do after Phase 1)

**Goal:** Replace legacy `params.x = value` declarations in `nextflow.config`
with a typed `params {}` block in `main.nf`.

### Current state (nextflow.config)

```groovy
params {
    trinity_fasta       = null
    samplesheet         = null
    species_label       = 'species_X'
    sortmerna_db_urls   = [...]
    sortmerna_db_dir    = null
    mmseqs2_swissprot   = null
    // ...
    filter_taxon        = 33090
    td2_min_orf_length  = 90
    td2_strand_specific = true
    max_cpus            = 16
    max_memory          = 128.GB
    max_time            = 168.h
    outdir              = './results'
}
```

### Typed equivalent (in main.nf, before `workflow {}`)

```groovy
params {
    // --- Input ---
    trinity_fasta:       Path                // required (no default)
    samplesheet:         Path                // required
    species_label:       String  = 'species_X'

    // --- SortMeRNA ---
    sortmerna_db_urls:   List<String> = [
        'https://raw.githubusercontent.com/.../rfam-5.8s-database-id98.fasta',
        // ...
    ]
    sortmerna_db_dir:    Path?               // optional, nullable

    // --- Databases ---
    mmseqs2_swissprot:   Path                // required
    mmseqs2_pfam:        Path                // required
    mmseqs2_eggnog:      Path                // required
    mmseqs2_taxonomy_db: Path                // required
    eggnog_annotations:  Path                // required
    busco_lineage:       String              // required
    diamond_db:          Path                // required

    // --- Cluster ---
    unix_group:          String? = null
    orion_exclude_nodes: String? = null

    // --- Search params ---
    mmseqs2_search_sens: Float   = 7.0
    filter_taxon:        Integer = 33090     // Viridiplantae

    // --- TD2 ---
    td2_min_orf_length:  Integer = 90
    td2_strand_specific: Boolean = true

    // --- Resource caps ---
    max_cpus:            Integer    = 16
    max_memory:          MemoryUnit = 128.GB
    max_time:            Duration   = 168.h

    // --- Output ---
    outdir:              Path = './results'
}
```

### What this gives us

- **Required params validated at launch** — no more manual `if (!params.x) error` blocks
- **Type coercion from CLI** — `--filter_taxon 33090` is auto-parsed as Integer
- **Self-documenting** — types serve as inline documentation
- **Language server validation** — type mismatches caught in editor before run

### Key decisions

| Param | Type | Rationale |
|-------|------|-----------|
| `sortmerna_db_dir` | `Path?` | Nullable — `null` means "download from URLs" |
| `unix_group`, `orion_exclude_nodes` | `String?` | Sometimes not set (non-Orion sites) |
| `filter_taxon` | `Integer` | NCBI taxon ID |
| `max_memory` | `MemoryUnit` | First-class unit type |
| `max_time` | `Duration` | First-class duration type |
| `outdir` | `Path` | Directory path |
| Container params | *Keep in config* | Not pipeline-level; stay in `nextflow.config` |

### Migration notes

- Typed params belong in `main.nf` (the script), not `nextflow.config`
- Config-only params (container paths, process selectors) stay in config
- Default values for required params can live in a `-params-file` or in a `test` profile
- The manual `if (!params.trinity_fasta) error` block can be deleted — typed params without defaults are required by the engine

### Estimated effort: ~1 hour, minimal risk (types are additive).

---

## Phase 3 — Typed Processes (preview, defer)

**Goal:** Migrate all 17 processes from legacy `val`/`path`/`tuple` qualifiers
to typed annotations with `name: Type` syntax.

### ⚠️ Status: Preview feature — syntax may change

Typed processes require:
1. `NXF_SYNTAX_PARSER=v2` (Phase 1)
2. `nextflow.preview.types = true` in every script that uses them
3. Rewrite of every `input:` and `output:` section

### What changes look like

**Before (legacy, current):**

```groovy
process MMSEQS2_TAXONOMY {
    input:
    path fasta
    path db
    val  sensitivity
    val  filter_taxon
    val  species_label

    output:
    path "supertranscripts_filtered.fasta", emit: fasta
    path "taxRes_lca.tsv",                  emit: lca_tsv
    path "taxonomy_filter_stats.txt",       emit: stats
    path "taxonomy_breakdown.tsv",          emit: breakdown

    script:
    // ...
}
```

**After (typed):**

```groovy
process MMSEQS2_TAXONOMY {
    input:
    fasta:         Path
    db:            Path
    sensitivity:   Float
    filter_taxon:  Integer
    species_label: String

    output:
    fasta:     Path     = file("supertranscripts_filtered.fasta")
    lca_tsv:   Path     = file("taxRes_lca.tsv")
    stats:     Path     = file("taxonomy_filter_stats.txt")
    breakdown: Path     = file("taxonomy_breakdown.tsv")

    script:
    // ...unchanged...
}
```

### Process-by-process conversion table

| Process | Inputs (legacy → typed) | Outputs (legacy → typed) | Complexity |
|---------|------------------------|--------------------------|------------|
| **SORTMERNA_INDEX** | `path(fastas)` → `fastas: Set<Path>` | `path "sortmerna_idx/"` → `file("sortmerna_idx/")` | Low |
| **SORTMERNA** | `tuple val(id), path(r1), path(r2)` → `(id, r1, r2): Tuple<String,Path,Path>` + `fastas: Set<Path>` + `index: Path` | 2× tuple → named tuple outputs | Medium |
| **SALMON_INDEX** | `path(fasta)` → `fasta: Path` | `path "salmon_idx/"` → `file("salmon_idx/")` | Low |
| **SALMON_QUANT** | `tuple val(id), path(r1), path(r2)` + `path(index)` + `val(suffix)` | `path "${id}_${suffix}/"` → `file(...)` | Medium |
| **CORSET** | `path(quant_dirs)` + `val(sample_conditions)` + `val(species_label)` | 2× `path` → 2× `file(...)` | Medium (the `sample_conditions` list-of-maps type is `List<Map>`) |
| **LACE** | `path(fasta)` + `path(clusters)` + `val(species_label)` | `path "supertranscripts.fasta"` → `file(...)` | Low |
| **MMSEQS2_TAXONOMY** | 3× `path` + 2× `val` | 4× `path` → 4× `file(...)` | Medium |
| **DIAMOND_BLASTX** | `path(fasta)` + `path(db)` + `val(species_label)` | `path "*.tsv"` → `file(...)` | Low |
| **CORRECT_FRAMESHIFTS** | `path(fasta)` + `path(tsv)` + `val(species_label)` | 2× `path` → 2× `file(...)` | Low |
| **TD2_LONGORFS** | `path(fasta)` + `val(species_label)` | 4× `path` → 4× `file(...)` | Medium (dir output) |
| **MMSEQS2_SEARCH** (×2) | `path(query_pep)` + `path(db)` + `val(tag)` | `path "*.m8"` → `file(...)` | Low |
| **TD2_PREDICT** | 5 inputs | 5 outputs | Medium |
| **SELECT_BEST_ORF** | 4 inputs | 3 outputs | Low |
| **BUSCO** (×2) | `path(fasta)` + `val(label)` + `val(suffix)` | 2 outputs | Low |
| **TRANSANNOT** | 6 inputs | 1 output | Low |
| **VALIDATE_IDS** | 3 inputs | 1 output | Low |
| **THINNING_REPORT** | 13 inputs | 1 output | Medium (many inputs) |

### Special considerations

1. **`stageAs` for name patterns:** Legacy `path "name_${var}", stageAs: '...'`
   must become explicit `stage:` directives in typed processes.

2. **`stageInMode 'copy'`:** TD2_LONGORFS and TD2_PREDICT use `stageInMode 'copy'`.
   This stays as a directive, unchanged.

3. **Collect patterns in outputs:** Patterns like `path "*.m8"` become `file("*.m8")`.
   For multiple files, use `files("*.m8")` → `Set<Path>`.

4. **Directory outputs:** `path "busco_${suffix}/"` — the trailing `/` matters.
   In typed syntax: `file("busco_${suffix}/")`.

5. **Tuple outputs (SORTMERNA):** emit tuples `[sample_id, r1, r2]` — in typed
   syntax this becomes named tuple outputs with an explicit `Tuple<String,Path,Path>` type.

6. **`collect()` inputs:** Several processes receive collected channels (e.g.,
   `CORSET` gets `path(quant_dirs)` from `.collect()`). In typed syntax,
   `Set<Path>` absorbs collected path channels naturally.

### Why defer

- **Preview feature**: The `nextflow.preview.types` flag may introduce breaking
  syntax changes in 26.04. Rewriting 17 processes now risks needing a second
  rewrite.
- **High effort**: Every process input/output must be manually converted and
  tested.
- **Cache invalidation**: Changing process definitions invalidates the Nextflow
  task cache — all tasks would re-run, even with `-resume`.
- **No functional benefit yet**: The type checker is language-server only; it
  doesn't affect runtime behavior.

### When to adopt

- When `nextflow.preview.types` exits preview (expected ~26.04–26.10)
- When the VS Code language server + `nextflow lint` can auto-convert (already
  partially works via "Convert script to static types" command)
- When no active runs depend on the current task cache

---

## Recommended Execution Order

```
Phase 1 — Strict Syntax          ← DO NOW (30 min, zero risk)
  ├── 1.1  Move workflow.onComplete → onComplete: section
  ├── 1.2  Move extractCondition() to lib/Utils.groovy
  ├── 1.3  Set NXF_SYNTAX_PARSER=v2 in run scripts
  └── 1.4  Validate with nextflow lint + -preview

Phase 2 — Typed Params           ← DO NEXT (1 hour, low risk)
  ├── 2.1  Add typed params {} block to main.nf
  ├── 2.2  Remove legacy params from nextflow.config (keep container params)
  ├── 2.3  Delete manual validation block
  └── 2.4  Test with nextflow run -preview + actual run

Phase 3 — Typed Processes        ← DEFER until preview exits
  ├── 3.1  Set nextflow.preview.types = true
  ├── 3.2  Convert processes one module at a time
  ├── 3.3  Run nextflow lint after each conversion
  ├── 3.4  Full integration test (fresh run, no -resume)
  └── 3.5  Remove nextflow.preview.types when stable
```

---

## Risk Matrix

| Change | Risk | Cache impact | Rollback |
|--------|------|-------------|----------|
| Phase 1 (strict syntax) | Very low | None | Unset `NXF_SYNTAX_PARSER` |
| Phase 2 (typed params) | Low | None (params don't affect task hash) | Restore `nextflow.config` params |
| Phase 3 (typed processes) | **Medium** | **Full cache invalidation** | Revert module files, restore cache |

---

## References

- [Preparing for strict syntax](https://www.nextflow.io/docs/latest/strict-syntax.html)
- [Typed processes](https://www.nextflow.io/docs/latest/process-typed.html)
- [Standard library types](https://www.nextflow.io/docs/latest/reference/stdlib-types.html)
- [Migrating to static types (tutorial)](https://www.nextflow.io/docs/latest/tutorials/static-types.html)
- [Migrating to workflow outputs](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
