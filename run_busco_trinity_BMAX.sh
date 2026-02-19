#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=busco-trinity-BMAX
#SBATCH --output=%x_%A.out
#SBATCH --chdir=/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/BMAX

# ── Purpose: Run BUSCO on the raw Trinity assembly (transcriptome mode)
#    to establish baseline completeness before pipeline processing.

module load singularity

unset BASH_ENV

TRINITY_FASTA=/mnt/project/FjellheimLab/martpali/AnnualPerennial/assemblies/BMAX-Trinity1/BMAX-Trinity.fasta
LINEAGE=/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db/busco_downloads/lineages/poales_odb12
OUTDIR=/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/BMAX/busco_trinity

mkdir -p "$OUTDIR"
cd "$OUTDIR"
rm -rf busco_trinity_BMAX

sg fjellheimlab -c "
singularity exec \
    --bind /mnt/project,/mnt/users \
    docker://ezlabgva/busco:v6.0.0_cv1 \
    busco \
        -i $TRINITY_FASTA \
        -m transcriptome \
        -l $LINEAGE \
        -o busco_trinity_BMAX \
        -c 32 \
        --offline
"
