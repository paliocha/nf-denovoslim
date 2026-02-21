# Changelog

All notable changes to nf-denovoslim are documented in this file.

## [Unreleased]

### Changed
- **Taxonomy filter fixed and tightened**:
  - Fixed critical bug: `filtertaxdb` + `createsubdb` was passing ALL sequences
    through (filtertaxdb zeros data but keeps keys; added `awk '$3 != 1'`
    intermediate step per MMseqs2 wiki)
  - Switched from `35493` (Streptophyta) + `2759` (Eukaryota) to `33090`
    (Viridiplantae) + `0` (no-hit only) — now removes fungi, nematodes,
    oomycetes, bacteria, archaea, viruses (~27% of BMAX SuperTranscripts)
  - Added taxonomy lineage breakdown table to filter output and thinning report
- **Switched from UniRef90 to UniRef50** for both taxonomy classification and
  DIAMOND frameshift correction — ~3-4× fewer sequences (52M vs 188M) gives
  proportionally faster searches with negligible effect on kingdom-level
  classification or frameshift correction accuracy
- **DIAMOND blastx optimised** — expected ~4–6× speedup (5.8 h → ~1–1.5 h):
  - `--strand plus`: skip reverse-strand translation (correction script already
    discards reverse hits); halves translated query volume
  - `--iterate --sensitive`: fast pre-screen at default sensitivity, only
    unmatched queries proceed to full sensitive search
  - `-g 512`: cap gapped extensions to top 512 ungapped-score targets;
    trims redundant Smith-Waterman work beyond `--top 1` best hit
  - `-b 4 -c 1`: larger reference blocks (34 → ~17) and single index pass
    (4 → 1); reduces per-iteration overhead by ~8×
  - `--tmpdir .`: ensure temp files go to node-local SSD
  - Memory raised 32 → 48 GB to support larger block size
- **Lace 2.0.0** — upgraded from patched v1.14.1 to complete rewrite
  ([paliocha/Lace](https://github.com/paliocha/Lace)) with:
  - minimap2 replaces BLAT (no licence restrictions)
  - Block-level directed graph replaces coloured de Bruijn graph
  - ProcessPoolExecutor parallelism (replaces GNU Parallel)
  - Python 3.12+ (tested on 3.12, 3.13 and 3.14)
  - **211× speedup** on BMAX dataset (1.17 M transcripts):
    v1.14.1 5 h 32 min → v2.0.0 94.5 s; identical cluster count
- **Corset 1.10** — upgraded from stock Corset 1.09 to OpenMP-parallelised
  fork ([paliocha/Corset](https://github.com/paliocha/Corset)) with:
  - Adjacency-list merge (O(degree) instead of O(n) per merge step)
  - Parallel distance recomputation (`#pragma omp parallel for`)
  - `FlatDistMap` open-addressing hash map (replaces `std::unordered_map`)
  - Early-out per-sample distance computation
  - SSE2 block intersection in `get_dist()`
  - **270× speedup** on FPRA dataset (40 samples, 793 K transcripts):
    v1.09 11 h 11 min → v1.10 3 min 21 s; 95.3 % identical cluster membership
  - **55× speedup** on BMAX dataset (40 samples, 460 K transcripts):
    v1.09 2 h 34 min → v1.10 2 min 47 s; 95.5 % identical cluster membership
- CORSET process default raised to 32 CPUs / 128 GB, `OMP_NUM_THREADS`
  exported automatically from `task.cpus`
- Container renamed from `corset_omp` to `corset`
- Container SIF tracked via Git LFS

### Fixed
- README clarified that the full unfiltered Trinity assembly must go into
  Salmon to preserve multi-mapping signal for Corset
- `maxForks` adjustments in base.config

## [0.3.0] — 2025-12-01

### Added
- BUSCO Trinity assessment (transcriptome mode) running in parallel with
  the main pipeline for baseline QC
- TransAnnot 4.0.0 functional annotation (SwissProt + Pfam + eggNOG)
- Thinning report summarising transcript counts at each pipeline stage

### Changed
- Taxonomy filtering switched from TrEMBL to UniRef90 (`--filter_taxon 35493`)
- Lace runs on node-local SSD to avoid NFS I/O bottleneck from per-cluster
  FASTA files
- TD2 updated to v1.0.8

### Fixed
- Lace container patched for NetworkX 3.x (`.node` → `.nodes`) and
  matplotlib < 3.6 pinning
- Diamond container lacked Python 3 — split into Diamond + Python processes
- Channel join bug dropping 39/40 samples (`9e0f376`)
- `bc` command not found in MMSEQS2_TAXONOMY

## [0.2.0] — 2025-09-01

### Added
- Frameshift correction via Diamond blastx + Python post-processing
- MMseqs2 taxonomy filtering (keep Streptophyta)
- Chunked subworkflows for frameshift correction and MMseqs2 search/taxonomy

### Changed
- Resource allocations optimised: escalating memory retries, `--split-memory-limit`
- SORTMERNA parallelisation increased (`maxForks 10 → 15`)

### Fixed
- DSL2 channel-fork bug starving SALMON_QUANT_INITIAL
- Singleton process outputs converted to value channels
- Scratch disabled globally (`.command.sh not found` on NFS)
- CORSET scratch disabled (40 dirs overflow local disk)

## [0.1.0] — 2025-06-01

### Added
- Initial pipeline: SortMeRNA → Salmon → Corset → Lace → TD2 → MMseqs2 →
  Select Best ORF → Salmon (gene-level) → BUSCO
- DSL2 modular structure with per-process containers
- Orion HPC profile (`-profile apptainer,orion`)
- Samplesheet-driven input with automatic condition extraction
- tximport-compatible Salmon quantification output
