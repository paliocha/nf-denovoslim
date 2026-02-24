# nf-denovoslim Pipeline Refactoring: Lace Removal, MetaEuk Integration, OrthoFinder-Ready Output

## Context

The nf-denovoslim pipeline thins Trinity de novo transcriptomes into non-redundant gene sets.
The current pipeline uses Lace (unmaintained since ~2017) to build SuperTranscripts, and TD2
alone for ORF prediction. This produces suboptimal results:

- **FPRA**: Trinity BUSCO 95.5% → Final protein BUSCO 56.8% (38.7pp loss)
- **BMED**: Trinity BUSCO 94.2% → Final protein BUSCO 76.7% (17.5pp loss)
- MetaEuk test on FPRA finds 34,260 genes that TD2 misses entirely
- MetaEuk alone achieves 60.0% BUSCO completeness (vs TD2's 56.8%) with 5.6pp fewer missing

## Aims

1. **Maintain or improve BUSCO completeness** — final protein BUSCO must be ≥ current scores
2. **Decrease duplicated BUSCOs** — MetaEuk alone has 32.3% duplication; target ≤20%
3. **Produce OrthoFinder-ready output** — best/longest ORFs and protein files per species

## New Pipeline Architecture

```
Trinity ─► BUSCO_TRINITY (baseline, unchanged)
         ─► SortMeRNA ─► Salmon initial quant ─► CORSET
                                                    │
                                              SELECT_REP (longest transcript per cluster)
                                                    │
                                              MMSEQS2_CLUSTER (95% nt identity dedup)
                                                    │
                                              MMSEQS2_TAXONOMY (Viridiplantae filter)
                                                    │
                                              DIAMOND + CORRECT_FRAMESHIFTS
                                                    │
                              ┌──────────────────────┴──────────────────────┐
                              │                                             │
                        TD2 branch                                   MetaEuk branch
                     TD2_LONGORFS                                  METAEUK_PREDICT
                    ┌──────┴──────┐                                       │
              MMSEQS2_SEARCH   MMSEQS2_SEARCH                    (internal SwissProt
              (SwissProt)      (Pfam)                              homology search)
                    └──────┬──────┘                                       │
                     TD2_PREDICT                                 PSAURON_METAEUK
                      (+ PSAURON)                               (re-score MetaEuk
                           │                                     predictions)
                    SELECT_BEST_ORF                                       │
                      (TD2 best)                              SELECT_BEST_METAEUK
                           │                                   (best per gene)
                           └──────────────┬───────────────────────┘
                                          │
                                   MERGE_PREDICTIONS
                                (longest protein per gene,
                                 min PSAURON threshold)
                                          │
                              ┌────────────┼────────────┐
                              │            │            │
                         BUSCO_QC    TRANSANNOT    SALMON_FINAL
                        (proteins)                (on representatives)
                              │            │            │
                              └────────────┴────────────┘
                                          │
                                   THINNING_REPORT
```

## Changes

### New modules

#### 1. `modules/select_representative.nf` — SELECT_REP

Pick the longest Trinity transcript per Corset cluster. For 66% singletons this is identity.
Write a FASTA of representative transcripts + a mapping TSV.

- **Input**: Trinity FASTA, Corset clusters file
- **Output**: `representatives.fasta` (emit: fasta), `rep_map.tsv` (emit: map)
- **Script**: `bin/select_representative.py` — read clusters TSV, index Trinity FASTA
  (BioPython SeqIO), pick longest per cluster, write FASTA with cluster ID as header
- **Container**: `biopython_container` (already exists)
- **Resources**: 1 CPU, 8 GB, 30 min

#### 2. `modules/mmseqs2_cluster.nf` — MMSEQS2_CLUSTER

Sequence-level dedup at 95% nucleotide identity. Catches near-identical transcripts that
Corset missed (different read sets, assembly duplicates). Uses `mmseqs easy-cluster`.

- **Input**: representatives FASTA
- **Output**: `representatives_dedup.fasta` (emit: fasta), `cluster_stats.txt` (emit: stats)
- **Script**:
  ```bash
  mmseqs easy-cluster representatives.fasta clust tmp \
      --min-seq-id 0.95 -c 0.8 --cov-mode 0 --threads ${task.cpus}
  # Extract representative sequences from cluster results
  # Report: N input → N output, N collapsed
  ```
- **Container**: `mmseqs2_container` (already exists)
- **Resources**: 8 CPU, 32 GB, 2h

#### 3. `modules/metaeuk.nf` — METAEUK_PREDICT

Run MetaEuk easy-predict against SwissProt. Select best protein per gene
(highest alignment score → lowest evalue → longest).

- **Input**: corrected representatives FASTA, SwissProt DB path, species label
- **Output**: `metaeuk_raw.faa` (emit: faa), `metaeuk.gff` (emit: gff), `metaeuk_map.tsv` (emit: map)
- **Script**: `metaeuk easy-predict` + inline Python post-processing (from test_metaeuk_fpra.sh
  lines 71-122, already validated)
  ```bash
  metaeuk easy-predict ${fasta} ${swissprot_db} metaeuk_out tmp_metaeuk \
      --protein 1 --strand 2 -s 5.7 --min-length 30 \
      --metaeuk-eval 0.001 --metaeuk-tcov 0.3 --max-intron 50000 \
      --threads ${task.cpus}
  # Python: best protein per gene by score → evalue → length
  ```
- **Container**: `quay.io/biocontainers/metaeuk:7.bba0d80--pl5321h6a68c12_2` (biocontainer)
- **Resources**: 32 CPU, 64 GB, 12h
- **Note**: SwissProt DB is already in MMseqs2 format at
  `/mnt/project/glowberry/transannot/db/SwissProtDB`

#### 4. `modules/psauron_metaeuk.nf` — PSAURON_METAEUK

Run PSAURON scoring on MetaEuk protein predictions. This enables quality-based filtering
and fair comparison with TD2 predictions in the merge step.

- **Input**: MetaEuk protein FASTA, corrected representatives FASTA (nucleotide context for PSAURON)
- **Output**: `metaeuk_psauron.csv` (emit: scores)
- **Script**: PSAURON is bundled with TD2. Run it on MetaEuk proteins:
  ```bash
  # PSAURON needs nucleotide context — extract CDS regions from MetaEuk GFF
  # Then score the protein sequences
  psauron ${metaeuk_faa} > metaeuk_psauron.csv
  ```
  Need to verify PSAURON's exact CLI interface in the TD2 container.
- **Container**: `td2_container` (PSAURON is included in TD2 install)
- **Resources**: 1 CPU, 8 GB, 2h

#### 5. `modules/merge_predictions.nf` — MERGE_PREDICTIONS

Merge TD2 and MetaEuk predictions into a single best-protein-per-gene set.

- **Input**: TD2 best ORF FASTA + map, MetaEuk best ORF FASTA + map + PSAURON scores,
  species label
- **Output**: `{species}.faa` (emit: faa), `{species}_longest.faa` (emit: longest_faa),
  `merge_map.tsv` (emit: map), `merge_stats.txt` (emit: stats)
- **Script**: `bin/merge_predictions.py`
  ```
  For each gene (Corset cluster ID):
    1. Collect TD2 prediction (if any) with PSAURON score
    2. Collect MetaEuk prediction (if any) with PSAURON score
    3. Filter: discard predictions with PSAURON < 0.3 (configurable)
    4. If both pass filter: keep the LONGER protein (OrthoFinder priority)
    5. If only one passes: keep it
    6. If neither passes: gene has no protein (dropped)

  Output two protein files:
    - {species}.faa: all merged proteins (for BUSCO, annotation)
    - {species}_longest.faa: same content but guaranteed longest ORF per gene
      (alias, same file — the merge already picks longest)

  Output merge_map.tsv: gene_id, source (td2/metaeuk/both), protein_length,
    psauron_score, td2_length, metaeuk_length
  ```
- **Container**: `biopython_container`
- **Resources**: 1 CPU, 4 GB, 30 min

### Modified modules

#### 6. `main.nf` — Complete rewiring

- Remove LACE import and invocation
- Add imports: SELECT_REP, MMSEQS2_CLUSTER, METAEUK_PREDICT, PSAURON_METAEUK, MERGE_PREDICTIONS
- Wire MetaEuk branch in parallel with TD2 branch (both fed by CORRECT_FRAMESHIFTS.out.fasta)
- MERGE_PREDICTIONS feeds into BUSCO_QC, TRANSANNOT, VALIDATE_IDS
- SALMON_INDEX_FINAL + SALMON_QUANT_FINAL now index/quantify on CORRECT_FRAMESHIFTS.out.fasta
  (representative transcripts, not SuperTranscripts — same input, just renamed semantically)
- THINNING_REPORT updated inputs (no more supertranscripts FASTA, add merge stats)
- Add `params.metaeuk_container`, `params.min_psauron` to parameter block

#### 7. `nextflow.config` — New container + params

```groovy
metaeuk_container = 'quay.io/biocontainers/metaeuk:7.bba0d80--pl5321h6a68c12_2'

// Merge prediction params
min_psauron        = 0.3    // minimum PSAURON score for merged predictions
```

#### 8. `conf/base.config` — Resource allocations for new processes

Add process blocks for: SELECT_REP, MMSEQS2_CLUSTER, METAEUK_PREDICT, PSAURON_METAEUK,
MERGE_PREDICTIONS with resources specified above. Add publishDir for new outputs:
- `SELECT_REP` → `representatives/`
- `MERGE_PREDICTIONS` → `proteins/`

#### 9. `bin/thinning_report.py` — Updated statistics

- Remove SuperTranscript section (no more Lace)
- Add representative selection stats (N transcripts → N representatives)
- Add MMseqs2 dedup stats (N before → N after, % collapsed)
- Add MetaEuk vs TD2 merge stats:
  - Genes with TD2-only prediction
  - Genes with MetaEuk-only prediction
  - Genes with both (which was selected)
  - Mean/median protein length per source
  - Mean/median PSAURON per source
- Update BUSCO section to note method

#### 10. `modules/thinning_report.nf` — Updated inputs

Remove `supertranscripts_fasta` input, add `merge_stats` and `dedup_stats` inputs.

### Removed modules

#### 11. Remove Lace

- Remove `include { LACE }` from main.nf
- Remove `LACE(...)` invocation
- Keep `containers/lace/` directory but it's no longer used
- `lace_container` param can be removed from nextflow.config

### MetaEuk duplication handling (PSAURON re-scoring strategy)

The MetaEuk FPRA test shows 32.3% duplicated BUSCOs (vs TD2's 15.0%). Root causes:

1. **Fragmented predictions**: MetaEuk calls N-terminal and C-terminal halves separately
   when there's a gap in the assembly. Both hit the same BUSCO HMM → inflated duplication.
2. **True paralogs**: Polyploid Poales have many paralogs. MetaEuk recovers more complete
   proteins on paralog copies → more legitimately duplicated BUSCOs.

The PSAURON re-scoring strategy addresses cause 1 without discarding cause 2:

- **PSAURON filters fragments**: Short protein fragments from MetaEuk typically score low
  on PSAURON (which evaluates full-length coding potential). The min_psauron=0.3 threshold
  removes these while keeping genuine proteins.
- **Longest-per-gene merge reduces remaining duplication**: When both TD2 and MetaEuk predict
  for the same gene, picking the longer protein eliminates the shorter fragment.
- **True paralogs preserved**: Paralogs on different genes (different Corset clusters) are
  kept regardless — OrthoFinder needs them for gene tree inference.

Expected impact on FPRA BUSCO:
- Completeness: ≥60% (captures both TD2-only and MetaEuk-only BUSCOs)
- Duplication: ~15-20% (one prediction per gene, fragments filtered)
- Missing: ≤25% (union of both methods' recoveries)

### OrthoFinder output

The merged `{species}.faa` file is directly OrthoFinder-ready:
- One protein per gene (Corset cluster)
- Longest available ORF selected (from whichever predictor)
- Quality-filtered (PSAURON ≥ 0.3)
- Headers are gene IDs (Cluster-X.Y format), stable across runs

For multi-species OrthoFinder runs, each species produces one `.faa` file:
```
orthofinder -f /path/to/faa_dir/ -t 32
```

## Verification

1. **Unit test each new script**: Run on FPRA data (SuperTranscripts already exist)
   - `select_representative.py`: verify output FASTA has one sequence per cluster,
     headers match cluster IDs, lengths are max per cluster
   - `merge_predictions.py`: verify gene counts match union of TD2 + MetaEuk,
     no duplicate gene IDs, longer protein always selected

2. **FPRA dry run**: Execute the refactored pipeline on FPRA with `-resume`
   (reuses cached SortMeRNA, Salmon, Corset results). Compare:
   - Final protein count vs current (264,204)
   - BUSCO completeness vs current (56.8%)
   - BUSCO duplication vs current (15.0%)
   - Annotation rate vs current (34.9%)

3. **BUSCO regression check**: Trinity baseline BUSCO must be unchanged (95.5% FPRA,
   94.2% BMED). Final protein BUSCO must be ≥ current values.

4. **OrthoFinder compatibility**: Run `orthofinder -f` on the merged FAA files from
   all three species (FPRA, BMED, BMAX) to verify format acceptance.

## File summary

| Action | File | Description |
|--------|------|-------------|
| CREATE | `modules/select_representative.nf` | Longest transcript per cluster |
| CREATE | `modules/mmseqs2_cluster.nf` | 95% nt identity dedup |
| CREATE | `modules/metaeuk.nf` | MetaEuk easy-predict + best-per-gene |
| CREATE | `modules/psauron_metaeuk.nf` | PSAURON re-scoring of MetaEuk proteins |
| CREATE | `modules/merge_predictions.nf` | TD2 + MetaEuk merge (longest, min PSAURON) |
| CREATE | `bin/select_representative.py` | Pick longest transcript per cluster |
| CREATE | `bin/merge_predictions.py` | Merge logic with PSAURON filter |
| MODIFY | `main.nf` | Remove Lace, add new modules, rewire |
| MODIFY | `nextflow.config` | Add metaeuk_container, min_psauron param |
| MODIFY | `conf/base.config` | Resource allocations for new processes |
| MODIFY | `bin/thinning_report.py` | Update stats for new pipeline structure |
| MODIFY | `modules/thinning_report.nf` | Update input declarations |
| REMOVE | Lace usage in `main.nf` | Remove import + invocation (keep container) |
