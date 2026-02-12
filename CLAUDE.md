# Trinity Assembly Thinning Pipeline — Nextflow Implementation Plan

## Overview

A Nextflow DSL2 pipeline to collapse a fragmented Trinity *de novo* transcriptome assembly into a non-redundant gene set, producing:
1. **SuperTranscript FASTA** — one sequence per gene (for Salmon quantification)
2. **Protein FASTA** (`.faa`) — one best protein per gene (for OrthoFinder / TransAnnot)
3. **GFF3** — ORF coordinates on SuperTranscripts
4. **Mapping file** — `gene_id ↔ orf_id ↔ PSAURON score ↔ length`
5. Gene-level **Salmon quant.sf** with IDs matching the `.faa`

---

## Container Strategy

| Tool | Container source | Notes |
|------|-----------------|-------|
| **SortMeRNA** | `community.wave.seqera.io/library/sortmerna:4.3.7` | **Matches nf-core/rnaseq v3.22.2** |
| **Salmon** | `biocontainers/salmon:1.10.3--h6dccd9a_2` | **Matches nf-core/rnaseq v3.22.2** |
| **MMseqs2** | `quay.io/biocontainers/mmseqs2:<tag>` | Biocontainers |
| **Corset** | `quay.io/biocontainers/corset:1.09--h077b44d_6` | Biocontainers (bioconda) |
| **Lace** | `quay.io/biocontainers/lace:1.14.1--pyh5e36f6f_0` | Biocontainers (bioconda, includes BLAT) |
| **Diamond** | `quay.io/biocontainers/diamond:2.1.22--h13889ed_0` | Frameshift correction (blastx -F 15) |
| **TD2** | **Must build** (see below) | Not yet on biocontainers |
| **BUSCO** | `quay.io/biocontainers/busco:<tag>` | QC step |
| **TransAnnot** | `quay.io/biocontainers/transannot:<tag>` | Functional annotation (SwissProt + Pfam + eggNOG) |
| **AGAT** | [biocontainers](https://biocontainers.pro/tools/agat) | GFF/GTF annotation handling, filtering, QC |
| **Python/BioPython** | Custom or biocontainers | For the best-ORF selection script |

### Container Build Strategy (no Docker on HPC)

Docker is **not available** on this SLURM HPC. Custom containers use [`mambaorg/micromamba`](https://github.com/mamba-org/micromamba-docker) as the base image and are built into `.sif` files via Apptainer.

**Key `mambaorg/micromamba` patterns:**
- Base image: `mambaorg/micromamba:<version>` (e.g., `mambaorg/micromamba:2.5.0`). Python is NOT included by default.
- Install conda packages: `RUN micromamba install -y -n base -c conda-forge <packages> && micromamba clean --all --yes`
- Install pip packages: set `ARG MAMBA_DOCKERFILE_ACTIVATE=1` first, then use shell-form `RUN pip install ...`
- System packages (apt): switch to `USER root`, install, then `USER $MAMBA_USER`
- File copies: always use `--chown=$MAMBA_USER:$MAMBA_USER`
- **Apptainer usage:** `apptainer run` auto-activates the conda env; `apptainer exec` needs `_entrypoint.sh` prefix; `apptainer shell` needs `--shell /usr/local/bin/_apptainer_shell.sh`

### TD2 Container (must build)

TD2 is a pip package (`pip install TD2`) with a PSAURON neural net dependency. It also ships Perl scripts in `util/`. Create a Dockerfile:

```dockerfile
FROM mambaorg/micromamba:2.5.0

# Install Python + Perl via conda
RUN micromamba install -y -n base -c conda-forge \
       python=3.11 \
       perl \
    && micromamba clean --all --yes

# Enable conda env for subsequent RUN steps
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install TD2 from PyPI (includes PSAURON)
RUN pip install --no-cache-dir TD2

# Verify
RUN TD2.LongOrfs --help && TD2.Predict --help

LABEL maintainer="martin.paliocha@nmbu.no"
LABEL description="TD2 (TransDecoder 2) with PSAURON for ORF prediction"
```

Build `.sif` directly with Apptainer (no Docker needed):

```bash
# Build .sif from Dockerfile on HPC
apptainer build td2_1.0.8.sif docker-archive://<(docker save td2:1.0.8)
# OR build via a remote builder / local machine with Docker, then transfer:
# docker build -t td2:1.0.8 .
# docker save td2:1.0.8 -o td2_1.0.8.tar
# apptainer build td2_1.0.8.sif docker-archive://td2_1.0.8.tar
```

---

## Pipeline DAG (Process Dependency Graph)

```
                        (input reads)
                               │
                               ▼
                      SORTMERNA_INDEX
                 (build index from 8 rRNA DBs)
                               │
                               ▼
                         SORTMERNA
                    (per sample, remove rRNA)
                               │
                               ▼
                      filtered reads (non-rRNA)
                               │
      trinity_assembly.fasta   │
                │              │
                ┌──────────────┼───────────────┐
                ▼              ▼               ▼
         MMSEQS2_CLUSTER_NT   SALMON_INDEX  (filtered reads)
           (97% nt dedup)        │               │
                │                ▼               │
                │          SALMON_QUANT_INITIAL ◄┘
                │           (per sample, --hardFilter --dumpEq)
                │                │
                │                ▼
                └──────► CORSET ◄────────────┘
                     (hierarchical clustering on eq classes)
                               │
                               ▼
                           LACE
                     (one SuperTranscript per gene)
                               │
                               ▼
                      MMSEQS2_TAXONOMY
              (taxonomy + filtertaxdb → plant only)
                               │
                               ▼
                   FRAMESHIFT_CORRECTION
              (Diamond blastx -F 15 + fix frameshifts)
                               │
                ┌──────────────┼──────────────────┐
                ▼              ▼                  ▼
         TD2_LONGORFS    SALMON_INDEX_ST   (filtered reads)
                │              │                  │
                ▼              ▼                  │
         MMSEQS2_SEARCH  SALMON_QUANT_FINAL ◄────┘
         (vs SwissProt      (gene-level)
          + Pfam)
                │
                ▼
         TD2_PREDICT
         (--retain-mmseqs-hits)
                │
                ▼
         SELECT_BEST_ORF
         (one protein per gene,
          mapping file output)
                │
                ├──► species_X.faa
                ├──► best_orfs.gff3
                ├──► orf_to_gene_map.tsv
                │
        ┌───────┴───────────┐
        ▼                   ▼
   BUSCO_QC            TRANSANNOT
                  (.faa vs SwissProt
                   + Pfam + eggNOG7
                   in one step)
```

---

## Nextflow Process Specifications

### 0a. `SORTMERNA_INDEX` — Build rRNA index (once)

**Container:** `community.wave.seqera.io/library/sortmerna:4.3.7--b730cad73fc42b8e`
**Singularity:** `https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/919d9c8f5f2c3221a94efe96b81bde0c953c13ebb0a1eca6b690b90666006cad/data`

Matches nf-core/rnaseq v3.22.2 SortMeRNA module. Builds a reusable index from 8 rRNA reference databases (downloaded once from SortMeRNA v4.3.4 repo).

**Default rRNA databases** (same as nf-core/rnaseq):

| Database | Source |
|----------|--------|
| `rfam-5.8s-database-id98.fasta` | Rfam |
| `rfam-5s-database-id98.fasta` | Rfam |
| `silva-arc-16s-id95.fasta` | SILVA Archaea 16S |
| `silva-arc-23s-id98.fasta` | SILVA Archaea 23S |
| `silva-bac-16s-id90.fasta` | SILVA Bacteria 16S |
| `silva-bac-23s-id98.fasta` | SILVA Bacteria 23S |
| `silva-euk-18s-id95.fasta` | SILVA Eukaryota 18S |
| `silva-euk-28s-id98.fasta` | SILVA Eukaryota 28S |

URLs: `https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/<filename>`

```bash
sortmerna \
  --ref rfam-5.8s-database-id98.fasta \
  --ref rfam-5s-database-id98.fasta \
  --ref silva-arc-16s-id95.fasta \
  --ref silva-arc-23s-id98.fasta \
  --ref silva-bac-16s-id90.fasta \
  --ref silva-bac-23s-id98.fasta \
  --ref silva-euk-18s-id95.fasta \
  --ref silva-euk-28s-id98.fasta \
  --workdir . \
  --index 1
```

**Input:** 8 rRNA FASTA files
**Output:** SortMeRNA index directory (reused by all samples)

**Resources:** `cpus: 4, memory: 16.GB, time: 1.h`

---

### 0b. `SORTMERNA` — Remove rRNA reads (per sample)

**Container:** same as 0a

Runs **per sample** (parallelized by Nextflow). Uses the pre-built index from step 0a.

```bash
sortmerna \
  --ref rfam-5.8s-database-id98.fasta \
  --ref rfam-5s-database-id98.fasta \
  --ref silva-arc-16s-id95.fasta \
  --ref silva-arc-23s-id98.fasta \
  --ref silva-bac-16s-id90.fasta \
  --ref silva-bac-23s-id98.fasta \
  --ref silva-euk-18s-id95.fasta \
  --ref silva-euk-28s-id98.fasta \
  --reads ${reads_1} --reads ${reads_2} \
  --threads ${task.cpus} \
  --workdir . \
  --aligned rRNA_reads \
  --fastx \
  --other non_rRNA_reads \
  --paired_in \
  --out2 \
  --num_alignments 1 \
  -v \
  --index 0
```

**Flags (matching nf-core/rnaseq defaults):**
- `--num_alignments 1` : stop after first rRNA hit (fast)
- `--paired_in` : if one mate matches rRNA, both classified as rRNA
- `--out2` : write paired reads into separate R1/R2 output files
- `--index 0` : use pre-built index (from step 0a)
- `--fastx` : output in FASTQ format

**Input:** per-sample paired reads + SortMeRNA index
**Output:** `non_rRNA_reads_fwd.fq.gz`, `non_rRNA_reads_rev.fq.gz` (filtered reads for all downstream steps)

**Resources:** `cpus: 12, memory: 72.GB, time: 16.h` (nf-core `process_high`)

---

### 1. `MMSEQS2_CLUSTER_NT` — Remove near-identical transcripts

**Container:** `quay.io/biocontainers/mmseqs2:<tag>`

```bash
mmseqs easy-cluster ${trinity_fasta} trinity_clu_nt tmp_nt \
  --min-seq-id 0.97 -c 0.8 --cov-mode 1 --threads ${task.cpus}
```

**Input:** `trinity_assembly.fasta`
**Output:** `trinity_clu_nt_rep_seq.fasta` (representative sequences), `trinity_clu_nt_cluster.tsv`

**Resources:** `cpus: 8, memory: 16.GB, time: 2.h`

---

### 2. `SALMON_INDEX_INITIAL` — Index the deduplicated assembly

**Container:** `biocontainers/salmon:1.10.3--h6dccd9a_2` *(matches nf-core/rnaseq v3.22.2)*

```bash
salmon index -t ${deduped_fasta} -i salmon_idx --keepDuplicates -p ${task.cpus}
```

**Input:** `trinity_clu_nt_rep_seq.fasta`
**Output:** `salmon_idx/` directory

**Resources:** `cpus: 4, memory: 16.GB, time: 1.h`

---

### 3. `SALMON_QUANT_INITIAL` — Quantify each sample (with eq classes for Corset)

**Container:** `biocontainers/salmon:1.10.3--h6dccd9a_2` *(matches nf-core/rnaseq v3.22.2)*

Runs **per sample** (parallelized by Nextflow).

```bash
salmon quant -i ${salmon_idx} -l A \
  -1 ${reads_1} -2 ${reads_2} \
  --validateMappings \
  --seqBias --gcBias --posBias \
  --hardFilter \
  --dumpEq \
  -p ${task.cpus} \
  -o ${sample_id}_quant
```

**Input:** salmon index + per-sample read pairs
**Output:** `${sample_id}_quant/` directories (containing `aux_info/eq_classes.txt`, `quant.sf`)

**Critical:** Salmon must be run with `--hardFilter` (disables range factorization) and `--dumpEq` for Corset. Do NOT use `--rangeFactorizationBins` — Corset requires equivalence classes without range factorization.

**Resources:** `cpus: 8, memory: 16.GB, time: 4.h`

---

### 4. `CORSET` — Hierarchical transcript-to-gene clustering

**Container:** `quay.io/biocontainers/corset:1.09--h077b44d_6`

Corset uses Salmon equivalence classes and hierarchical clustering with condition-aware paralog splitting to group transcripts into genes.

```bash
corset \
  -i salmon_eq_classes \
  -g ${conditions} \
  -n ${sample_names} \
  -p corset \
  ${eq_classes_files}
```

**Input:** All `${sample_id}_quant/aux_info/eq_classes.txt` files + condition groupings
**Output:** `corset-clusters.txt` (transcript → gene cluster mapping) + `corset-counts.txt` (gene-level raw counts)

**Important notes:**
- The `-g` flag takes comma-separated condition labels (e.g., `T1_L,T1_L,T1_R,T2_L,...`)
- The `-n` flag takes comma-separated sample names matching the eq_classes.txt file order
- The `corset-clusters.txt` file is a two-column TSV: `transcript_id\tcluster_id`
- Cluster IDs are hierarchical (e.g., `Cluster-1234.0`); sub-clusters arise from condition-aware paralog splitting

**Resources:** `cpus: 4, memory: 64.GB, time: 8.h`

---

### 5. `LACE` — Build SuperTranscripts from Corset clusters

**Container:** `quay.io/biocontainers/lace:1.14.1--pyh5e36f6f_0` (includes BLAT dependency)

Lace constructs SuperTranscripts by aligning transcripts within each cluster using BLAT and merging them into a non-redundant consensus sequence.

```bash
Lace ${deduped_fasta} ${corset_clusters} \
  -t --cores ${task.cpus} -o lace_out

cp lace_out/SuperDuper.fasta supertranscripts.fasta
```

**Input:** deduplicated Trinity FASTA + `corset-clusters.txt`
**Output:** `supertranscripts.fasta` (one SuperTranscript per gene)

**Note:** Corset's clusters.txt format (`transcript\tcluster`) is natively compatible with Lace — no column reordering needed (unlike the former Grouper→Trinity conversion).

**Resources:** `cpus: 8, memory: 64.GB, time: 8.h`

---

### 5b. `MMSEQS2_TAXONOMY` — Classify and filter SuperTranscripts by taxonomy

**Container:** `quay.io/biocontainers/mmseqs2:<tag>` (matched by `MMSEQS2_.*` wildcard)

Uses MMseqs2's native taxonomy workflow with `filtertaxdb` for hierarchical NCBI taxonomy filtering. More robust than grepping lineage strings — `filtertaxdb` uses the full NCBI taxonomy tree to include all descendants of the target taxon.

```bash
# Create query DB + run taxonomy assignment
mmseqs createdb ${supertranscripts_fasta} queryDB
mmseqs taxonomy queryDB ${taxonomy_db} taxResult tmp_tax \
    --tax-lineage 1 -s 7.0 --threads ${task.cpus}

# Export LCA report + Krona-compatible report
mmseqs createtsv queryDB taxResult taxRes_lca.tsv
mmseqs taxonomyreport ${taxonomy_db} taxResult taxRes_report

# Filter by NCBI taxon ID (keeps all descendants)
mmseqs filtertaxdb ${taxonomy_db} taxResult filteredTaxResult \
    --taxon-list ${params.filter_taxon}

# Extract matching sequences
mmseqs createsubdb filteredTaxResult queryDB filteredDB
mmseqs createsubdb filteredTaxResult queryDB_h filteredDB_h
mmseqs convert2fasta filteredDB supertranscripts_filtered.fasta
```

**Input:** `supertranscripts.fasta` (Lace SuperTranscripts)
**Output:** `supertranscripts_filtered.fasta` (plant-only), `taxRes_lca.tsv`, `taxRes_report`, `taxonomy_filter_stats.txt`

**Key parameter:** `--filter_taxon 35493` (Streptophyta). Uses NCBI taxonomy hierarchy — all descendants are automatically included.

**Resources:** `cpus: 16, memory: 64.GB, time: 8.h`

---

### 5c. `FRAMESHIFT_CORRECTION` — Detect and correct assembly frameshifts

**Container:** `quay.io/biocontainers/diamond:2.1.22--h13889ed_0`

**Method:** Diamond blastx with frameshift-tolerant alignment (`-F 15`) + BTOP string parsing (based on Leder et al. 2021, *J Evol Biol* 34:138). Assembly-induced frameshifts (indels during Trinity assembly) are detected via Diamond's frameshift-aware protein alignment and corrected in-place:
- `\` (+1 frameshift, extra base inserted) → **delete 1 base**
- `/` (-1 frameshift, base missing) → **insert 1 N**

This prevents TD2 from predicting fragmented ORFs due to premature stop codons or frame breaks.

```bash
# 1. Run Diamond blastx with frameshift-tolerant alignment
diamond blastx \
    -F 15 \
    --sensitive \
    --top 1 \
    --min-score 50 \
    -d ${diamond_db} \
    -q ${supertranscripts_fasta} \
    --outfmt 6 qseqid qstart qend qlen qframe btop \
    -p ${task.cpus} \
    -o diamond_fs.tsv

# 2. Parse BTOP strings and correct frameshifts in-place
correct_frameshifts.py \
    ${supertranscripts_fasta} \
    diamond_fs.tsv \
    supertranscripts_corrected.fasta
```

**Input:** `supertranscripts_filtered.fasta` (plant-only SuperTranscripts)
**Output:** `supertranscripts_corrected.fasta`, `frameshift_stats.txt`
**Database:** `${params.diamond_db}` (SwissProt `.dmnd` file or protein FASTA — required parameter)

**Resources:** `cpus: 16, memory: 32.GB, time: 4.h`

**Notes:**
- Only forward-strand hits (frames 1, 2, 3) are corrected — SuperTranscripts from Lace are correctly oriented
- Sequences without Diamond hits pass through unchanged (novel, non-coding, or too divergent)
- Corrected sequences are used for both TD2 ORF prediction AND Salmon final quantification

---

### 6. `TD2_LONGORFS` — Extract candidate ORFs from SuperTranscripts

**Container:** `td2:1.0.8` (custom, see above)

```bash
TD2.LongOrfs -t ${supertranscripts_fasta} -S
```

**Flags:**
- `-S` : sense-strand only (SuperTranscripts have known orientation)
- Default minimum ORF: 90 aa

**Input:** `supertranscripts.fasta` (Lace SuperTranscripts)
**Output:** `supertranscripts.fasta.TD2_dir/longest_orfs.pep` (+ `.cds`, `.gff3`)

**Resources:** `cpus: 2, memory: 8.GB, time: 1.h`

---

### 7. `MMSEQS2_SEARCH_SWISSPROT` — Homology search vs SwissProt

**Container:** `quay.io/biocontainers/mmseqs2:<tag>`

```bash
# Database should be pre-downloaded and provided as a pipeline param
mmseqs easy-search \
  ${longest_orfs_pep} \
  ${swissprot_db} \
  swissprot_alnRes.m8 tmp_sp \
  -s 7.0 --threads ${task.cpus}
```

**Input:** `longest_orfs.pep` + SwissProt MMseqs2 database
**Output:** `swissprot_alnRes.m8`

**Resources:** `cpus: 16, memory: 32.GB, time: 4.h`

---

### 8. `MMSEQS2_SEARCH_PFAM` — Profile search vs Pfam

**Container:** `quay.io/biocontainers/mmseqs2:<tag>`

```bash
mmseqs easy-search \
  ${longest_orfs_pep} \
  ${pfam_db} \
  pfam_alnRes.m8 tmp_pfam \
  -s 7.0 --threads ${task.cpus}
```

**Input:** `longest_orfs.pep` + Pfam-A.full MMseqs2 profile database
**Output:** `pfam_alnRes.m8`

Note: Steps 7 and 8 run **in parallel** (no dependency between them).

**Resources:** `cpus: 16, memory: 32.GB, time: 4.h`

---

### 9. `TD2_PREDICT` — Final ORF prediction with homology support

**Container:** `td2:1.0.8` (custom)

```bash
# Concatenate homology hits
cat ${swissprot_m8} ${pfam_m8} > combined_alnRes.m8

TD2.Predict -t ${supertranscripts_fasta} \
  --retain-mmseqs-hits combined_alnRes.m8
```

**Input:** `trinity_genes.fasta` + `combined_alnRes.m8`
**Output:**
- `trinity_genes.fasta.TD2.pep`
- `trinity_genes.fasta.TD2.cds`
- `trinity_genes.fasta.TD2.gff3`
- `trinity_genes.fasta.TD2.bed`
- `trinity_genes.fasta.TD2_dir/psauron_score.csv`

**Resources:** `cpus: 4, memory: 16.GB, time: 2.h` (GPU optional, speeds up PSAURON)

---

### 10. `SELECT_BEST_ORF` — One protein per gene + mapping file

**Container:** Python with BioPython (custom lightweight image, or add to TD2 image)

This is a custom Python script. Key logic:

```python
#!/usr/bin/env python3
"""Select the single best ORF per SuperTranscript gene.

Outputs:
  - species_X.faa         : proteins with gene-level IDs (for TransAnnot + OrthoFinder + Salmon join)
  - best_orfs.gff3        : filtered GFF3 for retained ORFs only
  - orf_to_gene_map.tsv   : gene_id → orf_id → psauron_score → orf_length
"""
import csv, sys
from collections import defaultdict
from Bio import SeqIO

# --- Parse PSAURON scores ---
scores = {}
with open(sys.argv[1]) as f:  # psauron_score.csv
    reader = csv.reader(f)
    next(reader)  # skip header
    for row in reader:
        scores[row[0]] = float(row[1])

# --- Group ORFs by parent SuperTranscript ---
orfs_by_gene = defaultdict(list)
for rec in SeqIO.parse(sys.argv[2], "fasta"):  # .TD2.pep
    # TD2 headers: >GENE_ID.pN GENE_ID|m.N type:... len:...
    orf_id = rec.id
    # Parent gene = everything before the last .pN
    gene_id = orf_id.rsplit(".", 1)[0]
    psauron = scores.get(orf_id, 0.0)
    orfs_by_gene[gene_id].append((psauron, len(rec.seq), orf_id, rec))

# --- Select best: highest PSAURON, then longest ---
kept_orfs = set()
with open("orf_to_gene_map.tsv", "w") as mapf:
    mapf.write("gene_id\torf_id\tpsauron_score\torf_length\n")
    with open("species_X.faa", "w") as out:
        for gene_id in sorted(orfs_by_gene):
            orfs = orfs_by_gene[gene_id]
            orfs.sort(key=lambda x: (x[0], x[1]), reverse=True)
            best_psauron, best_len, best_orf_id, best_rec = orfs[0]
            # Write protein with GENE ID as header (not orf ID)
            out.write(f">{gene_id}\n{str(best_rec.seq)}\n")
            mapf.write(f"{gene_id}\t{best_orf_id}\t{best_psauron}\t{best_len}\n")
            kept_orfs.add(best_orf_id)

# --- Filter GFF3 ---
with open(sys.argv[3]) as gff_in, open("best_orfs.gff3", "w") as gff_out:
    for line in gff_in:
        if line.startswith("#"):
            gff_out.write(line)
            continue
        # Keep lines where any kept orf_id appears in the attributes
        if any(orf_id in line for orf_id in kept_orfs):
            gff_out.write(line)

n_genes = len(orfs_by_gene)
n_total = sum(len(v) for v in orfs_by_gene.values())
print(f"Selected {n_genes} best proteins from {n_total} total ORFs "
      f"({n_total - n_genes} secondary ORFs removed)")
```

**Input:**
- `trinity_genes.fasta.TD2_dir/psauron_score.csv`
- `trinity_genes.fasta.TD2.pep`
- `trinity_genes.fasta.TD2.gff3`

**Output:**
- `species_X.faa` — headers = gene IDs = SuperTranscript IDs = Salmon quant.sf Name column
- `best_orfs.gff3` — ORF coords on SuperTranscripts (only kept ORFs)
- `orf_to_gene_map.tsv` — the mapping file

**Resources:** `cpus: 1, memory: 4.GB, time: 10.min`

---

### 11. `SALMON_INDEX_FINAL` + `SALMON_QUANT_FINAL` — Gene-level quantification

**Container:** `biocontainers/salmon:1.10.3--h6dccd9a_2` *(matches nf-core/rnaseq v3.22.2)*

```bash
salmon index -t ${supertranscripts_fasta} -i st_salmon_idx -p ${task.cpus}
salmon quant -i st_salmon_idx -l A \
  -1 ${reads_1} -2 ${reads_2} \
  --validateMappings \
  --seqBias --gcBias --posBias \
  --rangeFactorizationBins 4 \
  -p ${task.cpus} \
  -o ${sample_id}_st_quant
```

**Output:** `quant.sf` where `Name` column = SuperTranscript gene IDs

---

### 12. `VALIDATE_IDS` — Sanity check: .faa IDs ⊂ Salmon gene IDs

**Container:** any with bash/awk

```bash
# Extract Salmon gene IDs
cut -f1 ${quant_sf} | tail -n+2 | sort > ids_salmon.txt

# Extract .faa protein IDs
grep "^>" ${faa} | sed 's/^>//' | sort > ids_faa.txt

# Proteins must be a subset of Salmon genes
ORPHANS=$(comm -23 ids_faa.txt ids_salmon.txt | wc -l)
if [ "$ORPHANS" -gt 0 ]; then
    echo "ERROR: ${ORPHANS} protein IDs not found in Salmon quant!" >&2
    comm -23 ids_faa.txt ids_salmon.txt >&2
    exit 1
fi

# Report genes without ORFs
NO_ORF=$(comm -23 ids_salmon.txt ids_faa.txt | wc -l)
echo "Genes with protein: $(wc -l < ids_faa.txt)"
echo "Genes without ORF (non-coding/fragmented): ${NO_ORF}"
echo "Total SuperTranscript genes: $(wc -l < ids_salmon.txt)"
```

---

### 13. `BUSCO_QC` — Assess completeness of final protein set

**Container:** `quay.io/biocontainers/busco:<tag>`

```bash
busco -i ${faa} -l ${busco_lineage} -m proteins \
  -o busco_final -c ${task.cpus}
```

**Input:** `species_X.faa`
**Output:** BUSCO summary

**Resources:** `cpus: 8, memory: 16.GB, time: 2.h`

---

### 14. `TRANSANNOT` — Functional annotation of final proteins

**Container:** `quay.io/biocontainers/transannot:<tag>` (available via [Bioconda](https://bioconda.github.io/recipes/transannot/README.html))

[TransAnnot](https://github.com/soedinglab/transannot) (Zelenskaia et al. 2024) annotates proteins against SwissProt, Pfam, and eggNOG in a single command using MMseqs2 internally. ~333x faster than EnTAP, ~18x faster than eggNOG-mapper.

```bash
# Create MMseqs2 query database from the final .faa
transannot createquerydb ${faa} queryDB tmp_createdb

# Annotate against all 3 databases (searches run in parallel internally)
transannot annotate queryDB \
  ${pfam_db} \
  ${eggnog_db} \
  ${swissprot_db} \
  ${species_label}_transannot.tsv \
  tmp_annotate \
  --no-run-clust
```

**Flags:**
- `--no-run-clust` : skip internal linclust clustering (our .faa is already one protein per gene)
- Default E-value cutoff: < 1e-5
- SwissProt: best hit only per query; Pfam/eggNOG: all non-overlapping domain hits retained

**Input:** `species_X.faa` + pre-built MMseqs2 databases (SwissProtDB, PfamDB, eggNOG7DB)

**Output:** `${species_label}_transannot.tsv` — 10-column TSV:
| Column | Description |
|--------|-------------|
| queryID | Protein/gene ID |
| targetID | Hit ID in reference DB |
| qstart | Query start position |
| qend | Query end position |
| description | Functional description |
| E-value | Expectation value |
| sequenceIdentity | % identity |
| bitScore | Bit score |
| typeOfSearch | "sequence" or "profile" |
| nameOfDatabase | "SwissProt", "Pfam", or "eggNOG" |

A single query can have multiple rows (one per non-overlapping hit across all 3 databases).

**Resources:** `cpus: 16, memory: 32.GB, time: 2.h`

---

## Pipeline Parameters (`nextflow.config`)

```groovy
params {
    // --- Input ---
    trinity_fasta       = null          // Path to Trinity.fasta
    samplesheet         = null          // CSV: sample_id,condition,reads_1,reads_2
    species_label       = 'species_X'   // Used for output .faa naming

    // --- Databases (must be pre-built MMseqs2 databases) ---
    mmseqs2_swissprot   = '/mnt/project/glowberry/transannot/db/SwissProtDB'
    mmseqs2_pfam        = '/mnt/project/glowberry/transannot/db/PfamDB'
    mmseqs2_eggnog      = '/mnt/project/glowberry/transannot/db/eggNOG7DB'
    mmseqs2_taxonomy_db = '.../db/UniProtTrEMBLtaxdb'  // TrEMBL for taxonomy (broader coverage)
    busco_lineage       = 'poales_odb12'

    // --- Diamond DB for frameshift correction (.dmnd or protein FASTA) ---
    diamond_db          = null  // REQUIRED — build: diamond makedb --in uniprot_sprot.fasta.gz -d swissprot

    // --- Clustering params ---
    mmseqs2_nt_id       = 0.97         // Nucleotide dedup identity
    mmseqs2_nt_cov      = 0.8          // Nucleotide dedup coverage
    mmseqs2_search_sens = 7.0          // MMseqs2 search sensitivity for homology

    // --- Taxonomy filter (NCBI taxon ID — includes all descendants) ---
    filter_taxon      = 35493             // Streptophyta

    // --- TD2 params ---
    td2_min_orf_length  = 90           // Minimum ORF length (aa)
    td2_strand_specific = true         // -S flag for TD2.LongOrfs

    // --- Output ---
    outdir              = './results'
}

// --- Container definitions ---
// IMPORTANT: Salmon version pinned to match nf-core/rnaseq v3.22.2
// See: https://github.com/nf-core/modules/blob/master/modules/nf-core/salmon/quant/main.nf
def salmon_container = workflow.containerEngine == 'singularity'
    ? 'https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2'
    : 'biocontainers/salmon:1.10.3--h6dccd9a_2'

def sortmerna_container = workflow.containerEngine == 'singularity'
    ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/919d9c8f5f2c3221a94efe96b81bde0c953c13ebb0a1eca6b690b90666006cad/data'
    : 'community.wave.seqera.io/library/sortmerna:4.3.7--b730cad73fc42b8e'

process {
    withName: 'SORTMERNA.*'        { container = sortmerna_container }
    withName: 'MMSEQS2_.*'         { container = 'quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_3' }
    withName: 'SALMON_.*'          { container = salmon_container }
    withName: 'CORSET'             { container = 'quay.io/biocontainers/corset:1.09--h077b44d_6' }
    withName: 'LACE'                  { container = 'quay.io/biocontainers/lace:1.14.1--pyh5e36f6f_0' }
    withName: 'FRAMESHIFT_CORRECTION' { container = 'quay.io/biocontainers/diamond:2.1.22--h13889ed_0' }
    withName: 'TD2_.*'                { container = '/path/to/td2_1.0.8.sif' }
    withName: 'BUSCO_QC'           { container = 'ezlabgva/busco:v6.0.0_cv1' }
    withName: 'TRANSANNOT'         { container = 'quay.io/biocontainers/transannot:4.0.0--h4ac6f70_0' }
    withName: 'SELECT_BEST_ORF'    { container = 'quay.io/biocontainers/biopython:1.81' }
    withName: 'VALIDATE_IDS'       { container = 'ubuntu:22.04' }
}

// --- Singularity/Apptainer profile ---
profiles {
    apptainer {
        apptainer.enabled    = true
        apptainer.autoMounts = true
        apptainer.cacheDir   = "${projectDir}/singularity_cache"
        docker.enabled       = false
    }
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
        singularity.cacheDir   = "${projectDir}/singularity_cache"
        docker.enabled         = false
    }
}

// --- SLURM executor (typical HPC) ---
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue    = 'normal'
        process.clusterOptions = '--account=your_account'
    }
}
```

---

## Nextflow Module Structure

```
trinity-thinning/
├── main.nf
├── nextflow.config
├── modules/
│   ├── sortmerna.nf
│   ├── mmseqs2_cluster_nt.nf
│   ├── salmon_index.nf
│   ├── salmon_quant.nf
│   ├── corset.nf
│   ├── lace.nf
│   ├── mmseqs2_taxonomy.nf
│   ├── supertranscripts.nf
│   ├── td2_longorfs.nf
│   ├── mmseqs2_search.nf          # Reusable: parameterized by DB
│   ├── td2_predict.nf
│   ├── select_best_orf.nf
│   ├── validate_ids.nf
│   ├── busco.nf
│   └── transannot.nf
├── bin/
│   └── select_best_orf.py         # The Python script from step 10
├── containers/
│   └── td2/
│       └── Dockerfile
├── conf/
│   ├── base.config                # Default resource allocations
│   ├── slurm.config
│   └── test.config                # Small test dataset config
└── assets/
    └── samplesheet_example.csv
```

---

## Key Implementation Notes for Claude Code

1. **Corset condition flags**: The `CORSET` process builds `-g` (conditions) and `-n` (sample names) flags dynamically from the samplesheet metadata, similar to how Grouper's YAML config was generated.

2. **Salmon `--hardFilter` + `--dumpEq`**: These are REQUIRED for Corset. The `--hardFilter` flag disables range factorization, which is incompatible with Corset's clustering. The final Salmon quant (on SuperTranscripts) does NOT need them.

3. **Corset output → Lace input**: Corset's `clusters.txt` has format `transcript_id\tcluster_id`. Lace natively reads this format — no column reordering needed.

4. **TD2 working directory**: TD2 creates `${input_fasta}.TD2_dir/` in the directory where the input fasta lives. In Nextflow, use `stageInMode: 'copy'` or symlink carefully so TD2 can write alongside the input.

5. **PSAURON score CSV**: Located at `${input_fasta}.TD2_dir/psauron_score.csv`. This file is critical for the best-ORF selection — make sure it's captured as a process output.

6. **ID consistency chain**: The entire pipeline's correctness hinges on:
   ```
   SuperTranscript FASTA headers
       = Salmon quant.sf Name column
       = species_X.faa FASTA headers (after stripping .pN suffix)
       = TransAnnot queryID column
       = orf_to_gene_map.tsv gene_id column
   ```
   The `VALIDATE_IDS` process enforces this.

7. **MMseqs2 databases**: Should be **pre-built** and provided as pipeline params, not downloaded at runtime. Pre-download with:
   ```bash
   mmseqs databases UniProtKB/Swiss-Prot swissprot tmp
   mmseqs databases Pfam-A.full pfam tmp
   ```

7b. **Diamond database** (`diamond_db`): Required for frameshift correction. Build from SwissProt FASTA:
   ```bash
   wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
   diamond makedb --in uniprot_sprot.fasta.gz -d swissprot
   # → swissprot.dmnd — pass as --diamond_db /path/to/swissprot.dmnd
   ```
   Alternatively, pass any protein FASTA file (Diamond will auto-build the DB at runtime).

8. **Apptainer from Docker**: For containers from Docker Hub/quay.io, Nextflow + Apptainer will auto-pull and convert. For the custom TD2 container, pre-build the `.sif` file and reference it by path.

9. **Resource estimates** are for a typical plant transcriptome (~100K-300K Trinity transcripts). Adjust in `conf/base.config` based on your actual assembly size.

10. **Resume-ability**: All processes should be idempotent. Use `publishDir` with `mode: 'copy'` for final outputs only. Intermediate files stay in Nextflow's work directory for `-resume` support.

11. **Salmon version pinning**: The pipeline uses **Salmon 1.10.3** (`biocontainers/salmon:1.10.3--h6dccd9a_2`), which is the exact version used by the current nf-core/rnaseq pipeline (v3.22.2) and nf-core/modules master. This ensures index compatibility if you later want to compare quantification results between this pipeline and nf-core/rnaseq runs, and avoids subtle differences in quasi-mapping behavior between Salmon versions. For Singularity/Apptainer, the pre-built image is available at `https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2`.

---

## Input Data Paths

### Target Species (3 only)

Pipeline runs for **BMAX**, **BMED**, and **FPRA** (grass species, Poales).

| Species | Assembly dir | Trinity FASTA | Samples |
|---------|-------------|---------------|---------|
| BMAX | `$PROJECTS/FjellheimLab/martpali/AnnualPerennial/assemblies/BMAX-Trinity1/` | `BMAX-Trinity.fasta` | 40 |
| BMED | `$PROJECTS/FjellheimLab/martpali/AnnualPerennial/assemblies/BMED-Trinity1/` | `BMED-Trinity.fasta` | 40 |
| FPRA | `$PROJECTS/FjellheimLab/martpali/AnnualPerennial/assemblies/FPRA-Trinity1/` | `FPRA-Trinity.fasta` | 40 |

Note: All three use the "1" variant assembly directories.

### Paired-End Reads & Existing nf-core/rnaseq Samplesheets

Reads at `$PROJECTS/FjellheimLab/martpali/AnnualPerennial/rawdata/trimmed-paired/{SPECIES}/`.

Existing nf-core/rnaseq samplesheets (CSV: `sample,fastq_1,fastq_2,strandedness`):
- `$HOME/AnnualPerennial/rnaseq/BMAX.2.nfrna.csv`
- `$HOME/AnnualPerennial/rnaseq/BMED.2.nfrna.csv`
- `$HOME/AnnualPerennial/rnaseq/FPRA.2.nfrna.csv`

**Sample naming convention:** `{SPECIES}{IndivID}_{Timepoint}_{Tissue}`
- Example: `FPRA18_T2_L` = FPRA species, individual 18, timepoint T2, leaf
- Timepoints: T1, T2, T3, T4, T5
- Tissues: L (leaf), R (root)
- All samples are `unstranded`

**Condition for Corset** = `{Timepoint}_{Tissue}` (e.g., `T2_L`, `T5_R`) — 10 conditions per species (5 timepoints x 2 tissues).

### Databases (pre-built MMseqs2)

All at `/mnt/project/glowberry/transannot/db/`:

| Database | Path prefix |
|----------|------------|
| SwissProt | `/mnt/project/glowberry/transannot/db/SwissProtDB` |
| Pfam | `/mnt/project/glowberry/transannot/db/PfamDB` |
| eggNOG v7 | `/mnt/project/glowberry/transannot/db/eggNOG7DB` |
| eggNOG v7 profiles | `/mnt/project/glowberry/transannot/db/eggNOG7_profiles` |

### BUSCO

Use **BUSCO v6** (`ezlabgva/busco:v6.0.0_cv1`) with **OrthoDB v12** and **Poales lineage** (`poales_odb12`).

---

## AGAT — Gene Annotation Handling

[AGAT](https://github.com/NBISweden/AGAT) (Another GTF/GFF Analysis Toolkit) is used for GFF3/GTF annotation manipulation, filtering, and QC. Container available via [biocontainers](https://biocontainers.pro/tools/agat).

**Key capabilities for this pipeline:**
- **Sanitize GFF3**: `agat_convert_sp_gxf2gxf.pl` — parses any GFF/GTF variant, fixes missing parents/IDs/phases, outputs clean GFF3
- **Filter by ORF size**: `agat_sp_filter_by_ORF_size.pl`
- **Keep longest isoform**: `agat_sp_keep_longest_isoform.pl`
- **Extract sequences**: `agat_sp_extract_sequences.pl` — extract protein/CDS/transcript FASTA from GFF3 + reference
- **Annotation statistics**: `agat_sp_statistics.pl`
- **Manage IDs/attributes**: `agat_sp_manage_IDs.pl`, `agat_sp_manage_attributes.pl`
- **Add functional annotation**: `agat_sp_manage_functional_annotation.pl` — import InterProScan, BLAST results
- **Filter by attribute**: `agat_sp_filter_feature_by_attribute_value.pl`
- **Merge annotations**: `agat_sp_merge_annotations.pl`

**Note:** `_sp_` tools load the full annotation into memory (OMNISCIENT parser) and auto-correct structure. `_sq_` tools process line-by-line (faster, less memory, no auto-fix).

---

## Setup Environment

```bash
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"

# Load singularity & activate Nextflow environment
module load singularity/rpm
micromamba activate $HOME/micromamba/envs/Nextflow

# For GitHub CLI (gh):
# micromamba activate $HOME/micromamba/envs/base
```

---

## Final Outputs

| File | Description | IDs match? |
|------|-------------|------------|
| `results/supertranscripts/supertranscripts.fasta` | One SuperTranscript per gene (all) | Gene IDs |
| `results/taxonomy/supertranscripts_filtered.fasta` | Plant-only SuperTranscripts | Gene IDs |
| `results/taxonomy/taxRes_lca.tsv` | MMseqs2 LCA taxonomy assignments | Gene IDs |
| `results/salmon_final/{sample}_st_quant/quant.sf` | Gene-level counts | Gene IDs ✓ |
| `results/proteins/species_X.faa` | One protein per gene | Gene IDs ✓ |
| `results/annotation/best_orfs.gff3` | CDS on SuperTranscripts | Gene IDs as seqnames ✓ |
| `results/annotation/orf_to_gene_map.tsv` | gene↔orf↔score↔length | Links all IDs |
| `results/clustering/corset-clusters.txt` | Transcript → gene map | — |
| `results/clustering/corset-counts.txt` | Gene-level raw counts | Gene IDs |
| `results/qc/busco_final/` | BUSCO completeness | — |
| `results/transannot/${species}_transannot.tsv` | Functional annotation (SwissProt + Pfam + eggNOG) | Gene IDs ✓ |
