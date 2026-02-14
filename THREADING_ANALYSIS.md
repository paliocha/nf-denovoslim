# Threading Support Analysis for nf-denovoslim Tools

## Summary

Surveyed TD2, CORSET, MMseqs2, and Lace documentation to verify threading/CPU parameter support.

## Tool-by-Tool Analysis

### 1. TD2.LongOrfs ✅ **SUPPORTS THREADING** (needs fix)

**Current status**: Module allocates 4 CPUs but **does NOT pass threading parameter**

**Threading parameter**: `--threads` or `-@`
```bash
-@ THREADS, --threads THREADS
                      number of threads to use, default=16
```

**Current module code** (td2_longorfs.nf:24-27):
```bash
TD2.LongOrfs \
    -t ${supertranscripts_fasta} \
    ${strand_flag} \
    -m ${params.td2_min_orf_length}
```

**REQUIRED FIX**:
```bash
TD2.LongOrfs \
    -t ${supertranscripts_fasta} \
    ${strand_flag} \
    -m ${params.td2_min_orf_length} \
    --threads ${task.cpus}
```

**Recommendation**:
- **Keep 4 CPUs** or increase to 8 CPUs in base.config
- **Add `--threads ${task.cpus}`** to the module script
- This will enable actual parallel processing

---

### 2. TD2.Predict ❌ **NO THREADING SUPPORT**

**Threading parameter**: None found in help output

**Current status**: Module allocates 4 CPUs but tool is single-threaded

**Recommendation**:
- **Reduce to 1 CPU** in base.config
- Frees 3 CPUs for other jobs
- No performance loss (tool doesn't use multiple cores)

---

### 3. Corset ❌ **NO THREADING SUPPORT**

**Threading parameter**: None (checked full help output)

**Current status**: Module allocates 4 CPUs but tool is single-threaded

**Recommendation**:
- **Reduce to 1 CPU** in base.config
- Frees 3 CPUs for other jobs
- No performance loss (tool doesn't use multiple cores)

---

### 4. Lace ✅ **SUPPORTS THREADING** (already correct)

**Threading parameter**: `--cores`
```bash
--cores CORES     The number of cores you wish to run the job on
                  (default = 1)
```

**Current module code** (lace.nf:27-32):
```bash
Lace \
    ${deduped_fasta} \
    ${corset_clusters} \
    -t \
    --cores ${task.cpus} \
    -o lace_out
```

**Status**: ✅ **Already correctly implemented**

**Current allocation**: 8 CPUs (appropriate for BLAT-based alignment)

---

### 5. MMseqs2 ✅ **SUPPORTS THREADING** (already correct)

**Threading parameter**: `--threads`
```bash
--threads INT     Number of CPU-cores used (all by default) [16]
```

**Current module code** - All MMseqs2 modules use:
```bash
--threads ${task.cpus}
```

**Status**: ✅ **Already correctly implemented**

**Current allocations**:
- MMSEQS2_CLUSTER_NT: 16 CPUs ✓
- MMSEQS2_TAXONOMY: 64 CPUs ✓
- MMSEQS2_SEARCH: 32 CPUs ✓

---

## Required Changes to Modules

### 1. Fix TD2_LONGORFS (modules/td2_longorfs.nf)

**Line 24-27**, change:
```bash
TD2.LongOrfs \
    -t ${supertranscripts_fasta} \
    ${strand_flag} \
    -m ${params.td2_min_orf_length}
```

**To**:
```bash
TD2.LongOrfs \
    -t ${supertranscripts_fasta} \
    ${strand_flag} \
    -m ${params.td2_min_orf_length} \
    --threads ${task.cpus}
```

---

## Required Changes to Resource Config (conf/base.config)

### Current (incorrect):
```groovy
withName: 'TD2_LONGORFS' {
    cpus   = { check_max( 4, 'cpus' ) }
    memory = { check_max( 16.GB * task.attempt, 'memory' ) }
    time   = { check_max( 4.h * task.attempt, 'time' ) }
}
withName: 'TD2_PREDICT' {
    cpus   = { check_max( 4, 'cpus' ) }
    memory = { check_max( 16.GB * task.attempt, 'memory' ) }
    time   = { check_max( 4.h * task.attempt, 'time' ) }
}
withName: 'CORSET' {
    cpus    = { check_max( 4, 'cpus' ) }
    memory  = { check_max( 64.GB * task.attempt, 'memory' ) }
    time    = { check_max( 8.h * task.attempt, 'time' ) }
}
```

### Recommended (optimized):
```groovy
withName: 'TD2_LONGORFS' {
    cpus   = { check_max( 8, 'cpus' ) }          // Increased from 4, now uses threading
    memory = { check_max( 16.GB * task.attempt, 'memory' ) }
    time   = { check_max( 2.h * task.attempt, 'time' ) }  // Reduced (faster with 8 threads)
}
withName: 'TD2_PREDICT' {
    cpus   = { check_max( 1, 'cpus' ) }          // Reduced from 4 (single-threaded)
    memory = { check_max( 8.GB * task.attempt, 'memory' ) }  // Reduced from 16GB
    time   = { check_max( 2.h * task.attempt, 'time' ) }     // Reduced from 4h
}
withName: 'CORSET' {
    cpus    = { check_max( 1, 'cpus' ) }         // Reduced from 4 (single-threaded)
    memory  = { check_max( 32.GB * task.attempt, 'memory' ) }  // Reduced from 64GB
    time    = { check_max( 2.h * task.attempt, 'time' ) }      // Reduced from 8h
}
```

---

## Impact Summary

| Tool | Before | After | Change | Impact |
|------|--------|-------|--------|--------|
| **TD2_LONGORFS** | 4 CPUs (unused) | 8 CPUs (used) | +4 CPUs, **ADD threading** | **2-4× faster** |
| **TD2_PREDICT** | 4 CPUs (wasted) | 1 CPU | -3 CPUs freed | No slowdown |
| **CORSET** | 4 CPUs (wasted) | 1 CPU | -3 CPUs freed | No slowdown |
| **LACE** | 8 CPUs (used) | 8 CPUs (used) | No change | Already optimal |
| **MMseqs2** | Various (used) | Various (used) | No change | Already optimal |

**Net result**:
- TD2_LONGORFS becomes 2-4× faster
- 2 CPUs freed overall (4 wasted freed, 4 more allocated to TD2_LONGORFS, -3-3+4 = -2)
- Actually, net is: freed 3+3 = 6 CPUs, used +4 for TD2 = **2 CPUs freed** plus **faster TD2_LONGORFS**

---

## Verification Commands

To verify threading is working after changes:

```bash
# During a run, check that TD2.LongOrfs is using multiple threads:
top -b -n 1 -u $USER | grep TD2

# Or check SLURM logs for CPU efficiency:
sacct -j <jobid> --format=JobID,JobName,AllocCPUs,Elapsed,TotalCPU,MaxRSS,State
```

Look for `TotalCPU` approaching `AllocCPUs × Elapsed` for multi-threaded jobs.

---

## Updated IMPROVE_PLAN.md Recommendations

The IMPROVE_PLAN.md should be updated with this corrected information:

**Phase 1: Quick Wins - Updated**

1. **TD2_LONGORFS**:
   - Increase to 8 CPUs (from 4)
   - **ADD `--threads ${task.cpus}` to module** ← Critical fix
   - Expected: 2-4× speedup from actual threading

2. **TD2_PREDICT**:
   - Reduce to 1 CPU (from 4) ← Correct, single-threaded
   - Frees 3 CPUs

3. **CORSET**:
   - Reduce to 1 CPU (from 4) ← Correct, single-threaded
   - Frees 3 CPUs

4. **LACE**:
   - Keep at 8 CPUs ← Already using `--cores` correctly
   - No changes needed

5. **MMseqs2**:
   - Keep current threading ← Already using `--threads` correctly
   - No changes needed
