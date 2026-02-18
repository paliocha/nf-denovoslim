#!/bin/bash
#
# run_species.sh — Launch nf-denovoslim for a species on NMBU Orion HPC
#
# Usage:
#   sbatch run_species.sh BMAX
#   sbatch run_species.sh BMED
#   sbatch run_species.sh FPRA
#
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --job-name=denovoslim
#SBATCH --output=%x_%A.out

set -euo pipefail

# ── Species from CLI argument ─────────────────────────────────────────
SPECIES="${1:?Usage: sbatch run_species.sh <SPECIES>  (e.g. BMAX, BMED, FPRA)}"

# ── Paths ─────────────────────────────────────────────────────────────
PIPELINE_DIR="$HOME/AnnualPerennial/nf-denovoslim"
LAUNCH_DIR="$PIPELINE_DIR/runs/${SPECIES}"
WORK_DIR="$PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/${SPECIES}/work"
TRINITY_FASTA="$PROJECTS/FjellheimLab/martpali/AnnualPerennial/assemblies/${SPECIES}-Trinity1/${SPECIES}-Trinity.fasta"
SAMPLESHEET="$PIPELINE_DIR/${SPECIES}.samplesheet.csv"
OUTDIR="$PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/${SPECIES}"

# ── Databases ─────────────────────────────────────────────────────────
DB_BASE="$PROJECTS/glowberry/transannot/db"
LOCAL_DB="$PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db"

# ── Environment ───────────────────────────────────────────────────────
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate Nextflow

module load Java
module load Anaconda3
module load singularity

# ── Launch ────────────────────────────────────────────────────────────
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

# Run under fjellheimlab group so all files count against the project quota
sg fjellheimlab -c "
nextflow run $PIPELINE_DIR/main.nf \
    -profile apptainer,orion \
    -resume \
    -ansi-log false \
    -w $WORK_DIR \
    --trinity_fasta $TRINITY_FASTA \
    --samplesheet $SAMPLESHEET \
    --species_label $SPECIES \
    --mmseqs2_swissprot $DB_BASE/SwissProtDB \
    --mmseqs2_pfam $DB_BASE/PfamDB \
    --mmseqs2_eggnog $DB_BASE/eggNOG7_profiles \
    --eggnog_annotations $DB_BASE/e7_as_e5_annotations.tsv \
    --mmseqs2_taxonomy_db $LOCAL_DB/UniRef90taxdb \
    --diamond_db $LOCAL_DB/uniref90.dmnd \
    --busco_lineage poales_odb12 \
    --filter_taxon 35493 \
    --unix_group fjellheimlab \
    --orion_exclude_nodes cn-37 \
    --outdir $OUTDIR
"
