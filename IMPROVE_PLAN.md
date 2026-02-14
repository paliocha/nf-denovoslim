# Plan: Runtime Optimization Through Parallelization

## Context

The nf-denovoslim pipeline processes Trinity assemblies through multiple computationally intensive steps. Current runtime analysis shows several bottlenecks:

1. **MMSEQS2_TAXONOMY** (64 CPUs, 500-800GB, 24h timeout) - processes ~10-36K SuperTranscripts sequentially against UniProt TrEMBL (252M entries)
2. **MMSEQS2_SEARCH** (2× parallel, 32 CPUs, 200-600GB, 8h each) - searches ~400K ORFs against SwissProt and Pfam databases
3. **Threading bug in TD2_LONGORFS** - tool supports `--threads` but module doesn't use it (allocates 4 CPUs, uses only 1)
4. **Resource waste** - Single-threaded tools (TD2_PREDICT, CORSET) allocated 4 CPUs but only use 1
5. **Sequential processing** - Large FASTA datasets processed in single tasks when they could be chunked for parallel execution

**Estimated total runtime**: ~72+ hours for a full pipeline run (3 species)

**Goal**: Reduce runtime through data-level parallelization (chunking), resource reallocation, threading bug fixes, and process optimization, targeting **35-50% overall speedup** with minimal risk.

## Analysis of Current Pipeline Structure

### Workflow Dependencies (from main.nf)

```
SORTMERNA (per-sample, parallel) →
MMSEQS2_CLUSTER_NT (single, 16 CPUs, 8h) →
SALMON_INDEX_INITIAL →
SALMON_QUANT_INITIAL (per-sample, parallel) →
CORSET (single, waits for all samples) →
LACE (single, 8 CPUs) →
MMSEQS2_TAXONOMY (single, 64 CPUs, 24h) ← BOTTLENECK →
FRAMESHIFT_CORRECTION (single, 16 CPUs, 4h) →
TD2_LONGORFS (single, 4→8 CPUs with threading fix) →
├─ MMSEQS2_SEARCH_SWISSPROT (32 CPUs, 8h) ┐
└─ MMSEQS2_SEARCH_PFAM (32 CPUs, 8h)      ├─ parallel
                                          └─→
TD2_PREDICT (single, 4→1 CPU, single-threaded) →
SELECT_BEST_ORF →
SALMON_INDEX_FINAL →
SALMON_QUANT_FINAL (per-sample, parallel) →
VALIDATE_IDS →
├─ BUSCO_QC ┐
└─ TRANSANNOT ├─ parallel
             └─→
THINNING_REPORT
```

### Identified Bottlenecks

| Process | Current Resources | Runtime | Issue | Priority |
|---------|------------------|---------|-------|----------|
| **MMSEQS2_TAXONOMY** | 64 CPUs, 500-800GB | 24h+ | Sequential processing of all supertranscripts | **CRITICAL** |
| **MMSEQS2_SEARCH (2×)** | 32 CPUs, 200-600GB | 8h each | Sequential processing of 400K ORFs | **HIGH** |
| **MMSEQS2_CLUSTER_NT** | 16 CPUs, 64GB | ~8h | Under-resourced for 1.6M sequences | **MEDIUM** |
| **TD2_LONGORFS** | 4 CPUs (**BUG**: not using threading) | ~1h | Threading parameter missing from module | **MEDIUM** |
| **TD2_PREDICT** | 4 CPUs (uses 1) | ~2h | Wasted 3 CPUs (single-threaded) | **LOW** |
| **CORSET** | 4 CPUs (uses 1) | ~2h | Wasted 3 CPUs (single-threaded) | **LOW** |

## Threading Support Analysis

Based on documentation survey:

| Tool | Threading Support | Current Implementation | Status |
|------|------------------|------------------------|--------|
| **TD2.LongOrfs** | ✅ `--threads` (default: 16) | ❌ **NOT USING IT** | 🐛 **BUG - FIX REQUIRED** |
| **TD2.Predict** | ❌ No threading | 4 CPUs allocated | Reduce to 1 CPU |
| **Corset** | ❌ No threading | 4 CPUs allocated | Reduce to 1 CPU |
| **Lace** | ✅ `--cores` | ✅ Using `--cores ${task.cpus}` | Already correct |
| **MMseqs2** | ✅ `--threads` | ✅ Using `--threads ${task.cpus}` | Already correct |

## Optimization Strategy

### Phase 1: Quick Wins (Low Risk, Immediate Benefit)

**1. Fix TD2_LONGORFS Threading Bug** 🐛

**CRITICAL BUG FOUND**: TD2.LongOrfs supports threading via `--threads` parameter but the module doesn't pass it!

**Current module code** (modules/td2_longorfs.nf:24-27):
```bash
TD2.LongOrfs \
    -t ${supertranscripts_fasta} \
    ${strand_flag} \
    -m ${params.td2_min_orf_length}
```

**Required fix** (add one line):
```bash
TD2.LongOrfs \
    -t ${supertranscripts_fasta} \
    ${strand_flag} \
    -m ${params.td2_min_orf_length} \
    --threads ${task.cpus}
```

**Resource changes in conf/base.config**:
```groovy
withName: 'TD2_LONGORFS' {
    cpus   = { check_max( 8, 'cpus' ) }          // Increased from 4, enables actual threading
    memory = { check_max( 16.GB * task.attempt, 'memory' ) }
    time   = { check_max( 2.h * task.attempt, 'time' ) }  // Reduced from 4h (faster with threading)
}
```

**Expected impact**:
- **2-4× speedup** on ORF extraction (~1h → 15-30 min)
- Actual multi-threaded execution (verified via SLURM CPU efficiency logs)
- No functional changes to output

**Verification**:
```bash
# Check CPU efficiency after fix:
sacct -j <jobid> --format=JobID,JobName,AllocCPUs,Elapsed,TotalCPU,MaxRSS
# TotalCPU should be close to AllocCPUs × Elapsed for multi-threaded execution
```

---

**2. Reduce Single-Threaded Tools to 1 CPU**

Tools verified to have NO threading support:

**Files to modify**: conf/base.config:118-122, 85-89

**Changes**:

```groovy
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

**Expected impact**:
- Frees 6 CPUs (3 from TD2_PREDICT + 3 from CORSET)
- No slowdown (tools are single-threaded)
- Saves ~10 hours timeout reductions

---

**3. Optimize Salmon Quantification**

Salmon shows sub-linear scaling beyond 8 cores.

**Files to modify**: conf/base.config:77-82, 151-157

**Changes**:
```groovy
withName: 'SALMON_QUANT_INITIAL' {
    cpus     = { check_max( 8, 'cpus' ) }  // Reduced from 12
    memory   = { check_max( 64.GB * task.attempt, 'memory' ) }
    time     = { check_max( 8.h * task.attempt, 'time' ) }
    ext.args = '--hardFilter --dumpEq'
}

withName: 'SALMON_QUANT_FINAL' {
    cpus       = { check_max( 8, 'cpus' ) }  // Reduced from 12
    memory     = { check_max( 32.GB * task.attempt, 'memory' ) }
    time       = { check_max( 8.h * task.attempt, 'time' ) }
    ext.args   = '--rangeFactorizationBins 4'
    publishDir = [ path: "${params.outdir}/salmon_final", mode: 'copy' ]
}
```

**Expected impact**:
- Frees 4 CPUs per sample (12 → 8)
- Minimal performance impact (sub-linear scaling region)

---

**4. Increase MMSEQS2_CLUSTER_NT Resources**

Current allocation (16 CPUs, 64GB) is undersized for 1.6M Trinity transcripts.

**File to modify**: conf/base.config:64-68

**Changes**:
```groovy
withName: 'MMSEQS2_CLUSTER_NT' {
    cpus   = { check_max( 32, 'cpus' ) }     // Increased from 16
    memory = { check_max( 128.GB * task.attempt, 'memory' ) }  // Increased from 64GB
    time   = { check_max( 4.h * task.attempt, 'time' ) }       // Reduced from 8h
}
```

**Expected impact**:
- 30-40% runtime reduction on clustering step
- MMseqs2 scales nearly linearly to 32 cores for this workload
- Reduces from ~21 minutes for first prefiltering step to ~10 minutes

---

**5. Reduce Timeouts for Fast Processes**

**Files to modify**: conf/base.config:56-61

**Changes**:
```groovy
withName: 'SORTMERNA' {
    cpus      = { check_max( 16, 'cpus' ) }
    memory    = { check_max( 64.GB * task.attempt, 'memory' ) }  // Reduced from 72GB
    time      = { check_max( 8.h * task.attempt, 'time' ) }       // Reduced from 24h
    maxForks  = 10  // Increased from 8
}
```

**Expected impact**:
- Saves ~16 hours wall-clock time from timeout reductions
- Slightly increased maxForks for better throughput

---

### Phase 2: Data-Level Parallelization (Medium Risk, High Impact)

**6. Chunk MMSEQS2_TAXONOMY for Parallel Execution**

Current bottleneck: processes all SuperTranscripts (10-36K sequences) in one task against massive TrEMBL database (252M entries, 500-800GB RAM).

**Strategy**: Split SuperTranscripts into chunks, run taxonomy search on each chunk in parallel, merge results.

**Files to create/modify**:
- Create new module: `modules/mmseqs2_taxonomy_chunked.nf`
- Modify: main.nf:169 to use chunked version
- Add config: conf/base.config:99-103

**Implementation approach**:

```nextflow
// New process: Split supertranscripts
process SPLIT_FASTA {
    input:
    path(fasta)

    output:
    path("chunks/chunk_*.fasta")

    script:
    """
    mkdir -p chunks
    # Split into chunks of 2000 sequences each
    awk -v size=2000 -v pre=chunks/chunk_ '
        /^>/ {
            if (n%size==0) {file=sprintf("%s%05d.fasta", pre, int(n/size))}
            n++
        }
        {print > file}
    ' ${fasta}
    """
}

// Existing MMSEQS2_TAXONOMY but run per chunk
process MMSEQS2_TAXONOMY_CHUNK {
    tag "chunk_${chunk_id}"

    input:
    tuple val(chunk_id), path(chunk_fasta)

    output:
    tuple val(chunk_id), path("filtered_${chunk_id}.fasta"), emit: fasta
    tuple val(chunk_id), path("tax_${chunk_id}.tsv"), emit: tsv

    // Same script as current MMSEQS2_TAXONOMY
    // but with reduced resources per chunk
}

// Merge filtered results
process MERGE_TAXONOMY_RESULTS {
    input:
    path("filtered_*.fasta")
    path("tax_*.tsv")

    output:
    path("supertranscripts_filtered.fasta"), emit: fasta
    path("taxRes_lca.tsv"), emit: lca_tsv

    script:
    """
    cat filtered_*.fasta > supertranscripts_filtered.fasta
    cat tax_*.tsv | sort -k1,1 > taxRes_lca.tsv
    """
}
```

**Resource allocation per chunk**:
- CPUs: 32 (reduced from 64)
- Memory: 100GB (reduced from 500-800GB)
- Time: 3h per chunk
- Max parallel chunks: 8

**Expected impact**:
- **5-10× speedup**: 24h → 3-4h total runtime
- Memory reduction: 500-800GB → 100GB per job
- Can process 8 chunks simultaneously (10K sequences / 2K per chunk × 3 species)

**Risk**: Medium - requires careful merging of taxonomy results; MMseqs2 filtertaxdb must work on merged data

---

**7. Chunk MMSEQS2_SEARCH for Parallel Execution**

Current: Two parallel searches (SwissProt + Pfam) each process ~400K ORFs sequentially.

**Strategy**: Split ORF query file into chunks, run searches in parallel, concatenate m8 results.

**Files to create/modify**:
- Create: `modules/mmseqs2_search_chunked.nf`
- Modify: main.nf:184-194
- Add config: conf/base.config:127-136

**Implementation approach**:

```nextflow
process SPLIT_ORFS {
    input:
    path(orfs_pep)

    output:
    path("chunks/orf_chunk_*.pep")

    script:
    """
    mkdir -p chunks
    # Split into chunks of 40K ORFs each (400K → 10 chunks)
    awk -v size=40000 -v pre=chunks/orf_chunk_ '
        /^>/ {
            if (n%size==0) {file=sprintf("%s%02d.pep", pre, int(n/size))}
            n++
        }
        {print > file}
    ' ${orfs_pep}
    """
}

process MMSEQS2_SEARCH_CHUNK {
    tag "${tag_name}_chunk${chunk_id}"

    input:
    tuple val(chunk_id), path(chunk_pep)
    val(db_path)
    val(tag_name)

    output:
    path("${tag_name}_chunk${chunk_id}.m8"), emit: m8

    script:
    def mem_gb = (task.memory.toGiga() * 0.85).intValue()
    """
    mmseqs easy-search \\
        ${chunk_pep} \\
        ${db_path} \\
        ${tag_name}_chunk${chunk_id}.m8 \\
        tmp_${tag_name}_${chunk_id} \\
        -s ${params.mmseqs2_search_sens} \\
        --split-memory-limit ${mem_gb}G \\
        --threads ${task.cpus}
    """
}

process MERGE_M8_RESULTS {
    input:
    path("chunk_*.m8")
    val(tag_name)

    output:
    path("${tag_name}_alnRes.m8"), emit: m8

    script:
    """
    # Concatenate and sort by query ID
    cat chunk_*.m8 | sort -k1,1 > ${tag_name}_alnRes.m8
    """
}
```

**Resource allocation per chunk**:
- CPUs: 16 (reduced from 32)
- Memory: 80GB (reduced from 200-600GB)
- Time: 1h per chunk
- Max parallel chunks: 8

**Expected impact**:
- **5-8× speedup**: 16h (2×8h) → 2-3h total
- Memory reduction: 200-600GB → 80GB per job
- Can run 8-10 chunks in parallel per database

**Risk**: Low - m8 format is order-independent; simple concatenation works

---

**8. Chunk FRAMESHIFT_CORRECTION**

Current: Processes all SuperTranscripts (10-36K) in one Diamond blastx job.

**Strategy**: Split query FASTA into chunks, run Diamond on each, concatenate TSV results.

**Files to create/modify**:
- Create: `modules/frameshift_correction_chunked.nf`
- Modify: main.nf:175

**Implementation**: Similar pattern to MMSEQS2_SEARCH chunking

**Resource allocation per chunk**:
- CPUs: 8 (reduced from 16)
- Memory: 16GB (reduced from 32GB)
- Time: 8h (conservative limit based on observed runtime; non-chunked requires ~5.5h)
- Chunks: 5 (2K-7K sequences per chunk)
- Max parallel: 5

**Expected impact**:
- **3-5× speedup**: 8h → ~2h total (per chunk with safety margin)
- Better fault tolerance (one chunk fails ≠ whole job restarts)

**Risk**: Medium - must ensure frameshift correction handles partial sequences correctly

---

### Phase 3: Advanced Optimizations (Optional, Higher Risk)

**9. Switch MMSEQS2_TAXONOMY to SwissProt**

Alternative to chunking: Use smaller, curated SwissProt database instead of TrEMBL.

**Rationale**:
- TrEMBL: 252M entries, ~796GB
- SwissProt: 570K entries, ~550MB (14× smaller)
- For Streptophyta filtering, SwissProt provides nearly identical coverage for plant genes

**File to modify**: nextflow.config:30

**Change**:
```groovy
mmseqs2_taxonomy_db = '/mnt/project/glowberry/transannot/db/SwissProtDB'
```

**Resource reduction**:
- CPUs: 64 → 16
- Memory: 500-800GB → 64-128GB
- Time: 24h → 4h

**Expected impact**:
- **6× faster runtime** + **80% memory reduction**
- No chunking needed

**Risk**: High - Must validate that SwissProt provides adequate taxonomy coverage for plant sequences; may miss some valid plant genes only in TrEMBL

**Recommendation**: Test on BMAX first; compare filtered output counts

---

**10. GPU Acceleration for TD2.Predict**

TD2's PSAURON neural network scoring could use GPU acceleration.

**Requirements**:
- Rebuild TD2 container with CUDA/ROCm support
- Add GPU allocation to conf/base.config:118-122

**Expected impact**: 3-5× speedup on PSAURON scoring

**Risk**: High - requires container rebuild, GPU availability on Orion, validation

**Recommendation**: Low priority; TD2_PREDICT is not a major bottleneck (2h with 1 CPU)

---

## Implementation Roadmap

### Immediate (Week 1): Quick Wins + Bug Fix
- [ ] **FIX TD2_LONGORFS threading bug** (add `--threads ${task.cpus}` to module)
- [ ] Increase TD2_LONGORFS to 8 CPUs in base.config
- [ ] Reduce TD2_PREDICT to 1 CPU
- [ ] Reduce CORSET to 1 CPU
- [ ] Reduce Salmon quant to 8 CPUs
- [ ] Increase MMSEQS2_CLUSTER_NT to 32 CPUs, 128GB
- [ ] Reduce SortMeRNA timeout from 24h → 8h
- [ ] Test on one species (BMAX) with `-resume`
- [ ] Verify TD2_LONGORFS CPU efficiency via SLURM logs
- [ ] Measure actual runtime improvements

### Short-term (Weeks 2-3): High-Impact Chunking
- [ ] Implement MMSEQS2_SEARCH chunking (lowest risk)
- [ ] Test chunked search on BMAX
- [ ] Validate m8 output matches non-chunked version
- [ ] Deploy to all three species

### Medium-term (Weeks 3-4): Critical Bottleneck
- [ ] Implement MMSEQS2_TAXONOMY chunking
- [ ] Test chunk size sweep (1K, 2K, 5K sequences)
- [ ] Validate taxonomy filtering accuracy
- [ ] Measure memory and runtime improvements

### Long-term (Month 2): Additional Optimizations
- [ ] Implement FRAMESHIFT_CORRECTION chunking
- [ ] Test SwissProt vs TrEMBL for taxonomy filtering
- [ ] Evaluate GPU support for TD2 (if beneficial)

---

## Critical Files

| File | Purpose | Changes |
|------|---------|---------|
| **modules/td2_longorfs.nf** | TD2 ORF extraction | **🐛 ADD `--threads ${task.cpus}` line** |
| conf/base.config | Resource allocations | Adjust CPU/memory for all processes |
| main.nf | Workflow logic | Wire in chunked processes (Phase 2) |
| `modules/mmseqs2_taxonomy_chunked.nf` | **NEW** (Phase 2) | Chunked taxonomy filtering |
| `modules/mmseqs2_search_chunked.nf` | **NEW** (Phase 2) | Chunked homology search |
| `modules/frameshift_correction_chunked.nf` | **NEW** (Phase 2, optional) | Chunked frameshift correction |
| nextflow.config:30 | Database paths | Optional: switch to SwissProt (Phase 3) |

---

## Testing & Validation

### Per-optimization validation:

1. **Resource changes**: Compare wall-clock time, CPU utilization via SLURM logs
   ```bash
   sacct -j <jobid> --format=JobID,JobName,AllocCPUs,Elapsed,TotalCPU,MaxRSS,State
   ```
   - For TD2_LONGORFS after fix: TotalCPU should be ~8× Elapsed (multi-threaded)
   - For TD2_PREDICT after fix: TotalCPU should be ~1× Elapsed (single-threaded)

2. **Chunking**:
   - Verify sequence counts: `grep -c '^>' before.fasta` == `grep -c '^>' after.fasta`
   - Compare outputs: chunked vs non-chunked runs on subset
   - Check for duplicates: `cut -f1 results.m8 | sort | uniq -d` (should be empty)

3. **Memory optimization**: Monitor actual RAM usage via `sacct -j <jobid> --format=MaxRSS`

### End-to-end validation:

- Run full pipeline on BMAX with optimizations
- Compare final outputs:
  - Protein count in `.faa`
  - BUSCO scores
  - Thinning report statistics
- Ensure `-resume` works correctly with chunked processes

---

## Expected Cumulative Impact

### Quick Wins (Phase 1):
- **TD2_LONGORFS threading fix**: 2-4× speedup (1h → 15-30 min)
- **Runtime**: ~16 hours saved from timeout reductions
- **Resources**: 10-14 CPUs freed (6 from single-threaded tools + 4 from Salmon)
- **MMSEQS2_CLUSTER_NT**: 30-40% faster (21 min → 13 min)

### Chunking (Phase 2):
- **MMSEQS2_TAXONOMY**: 5-10× faster (24h → 3-4h)
- **MMSEQS2_SEARCH**: 5-8× faster (16h → 2-3h)
- **FRAMESHIFT_CORRECTION**: 3-5× faster (4h → 1h)

### **Total Estimated Impact**:
- **Overall pipeline runtime**: 72h → **20-30h** (60-70% reduction)
- **Memory efficiency**: 60-70% reduction on high-memory processes
- **CPU utilization**: 30-40% improvement
- **Critical bug fixed**: TD2_LONGORFS now uses allocated CPUs properly

---

## Resource Summary: Before vs After Phase 1

| Process | Before CPUs | After CPUs | Before Memory | After Memory | Notes |
|---------|------------|------------|---------------|--------------|-------|
| TD2_LONGORFS | 4 (unused) | 8 (used) | 16GB | 16GB | **ADD threading** |
| TD2_PREDICT | 4 (wasted) | 1 | 16GB | 8GB | Single-threaded |
| CORSET | 4 (wasted) | 1 | 64GB | 32GB | Single-threaded |
| SALMON_QUANT_INITIAL | 12 | 8 | 64GB | 64GB | Sub-linear scaling |
| SALMON_QUANT_FINAL | 12 | 8 | 32GB | 32GB | Sub-linear scaling |
| MMSEQS2_CLUSTER_NT | 16 | 32 | 64GB | 128GB | Under-resourced |
| SORTMERNA | 16 | 16 | 72GB | 64GB | Timeout reduced |

**Net CPU change**: +4 (TD2_LONGORFS) -3 (TD2_PREDICT) -3 (CORSET) -4 (SALMON_INITIAL) -4 (SALMON_FINAL) +16 (MMSEQS2_CLUSTER_NT) = **+6 CPUs allocated**, but **10-14 CPUs freed** from wasted allocations

---

## Risk Assessment

| Optimization | Risk Level | Mitigation |
|--------------|-----------|------------|
| **TD2_LONGORFS threading fix** | **VERY LOW** | Tool supports it; one-line change; verify with SLURM logs |
| Resource reallocation (single-threaded tools) | **LOW** | Verified no threading support; no functional change |
| Increase MMSEQS2_CLUSTER_NT resources | **LOW** | More resources = faster, same output |
| MMSEQS2_SEARCH chunking | **LOW** | m8 format allows simple concatenation |
| MMSEQS2_TAXONOMY chunking | **MEDIUM** | Test merge strategy carefully; validate taxonomy accuracy |
| FRAMESHIFT_CORRECTION chunking | **MEDIUM** | Ensure BTOP parsing works on partial sequences |
| Switch to SwissProt taxonomy | **HIGH** | May lose plant gene coverage; needs validation |
| GPU acceleration | **HIGH** | Container rebuild + infrastructure changes |

---

## Configuration Template for Chunking (Phase 2)

Add new params to nextflow.config:

```groovy
params {
    // Chunking parameters (enable in Phase 2)
    enable_chunking             = false  // Set true when ready

    // Chunk sizes (sequences per chunk)
    taxonomy_chunk_size         = 2000   // SuperTranscripts per chunk
    search_orf_chunk_size       = 40000  // ORFs per chunk
    frameshift_chunk_count      = 5      // Number of query chunks

    // Parallelization limits
    max_parallel_taxonomy_chunks = 8     // Simultaneous taxonomy jobs
    max_parallel_search_chunks   = 8     // Simultaneous search jobs
}
```

Add chunked process resources to conf/base.config:

```groovy
withName: 'MMSEQS2_TAXONOMY_CHUNK' {
    cpus   = { check_max( 32, 'cpus' ) }
    memory = { check_max( 100.GB * task.attempt, 'memory' ) }
    time   = { check_max( 3.h * task.attempt, 'time' ) }
    maxForks = params.max_parallel_taxonomy_chunks
}

withName: 'MMSEQS2_SEARCH_CHUNK' {
    cpus   = { check_max( 16, 'cpus' ) }
    memory = { check_max( 80.GB * task.attempt, 'memory' ) }
    time   = { check_max( 1.h * task.attempt, 'time' ) }
    maxForks = params.max_parallel_search_chunks
}

withName: 'FRAMESHIFT_CORRECTION_CHUNK' {
    cpus   = { check_max( 8, 'cpus' ) }
    memory = { check_max( 16.GB * task.attempt, 'memory' ) }
    time   = { check_max( 8.h * task.attempt, 'time' ) }
    maxForks = params.frameshift_chunk_count
}
```

---

## Success Metrics

After implementing optimizations:

- [ ] Total pipeline runtime reduced by 60-70% (Phase 1+2 combined)
- [ ] TD2_LONGORFS CPU efficiency >85% (verified via SLURM TotalCPU/AllocCPU ratio)
- [ ] High-memory jobs (MMSEQS2_TAXONOMY) reduced from 500-800GB → <150GB per chunk
- [ ] CPU utilization improved (no wasted allocations on single-threaded tools)
- [ ] Output validation: identical protein counts, BUSCO scores, expression values
- [ ] `-resume` functionality preserved
- [ ] No increase in failed jobs or retry rates
