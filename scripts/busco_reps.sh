#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=busco_reps
#SBATCH --output=%x_%A.out
#SBATCH --chdir=/net/fs-2/scale/OrionStore/Home/martpali/AnnualPerennial/nf-denovoslim/runs/BMED

CONTAINER=/mnt/users/martpali/AnnualPerennial/nf-denovoslim/singularity_cache/ezlabgva-busco-v6.0.0_cv1.img
FASTA=/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/BMED/representatives/representatives.fasta
LINEAGE=poales_odb12
CPUS=16

module load singularity

apptainer exec \
    -B /mnt/project:/mnt/project \
    -B /net/fs-2/scale:/net/fs-2/scale \
    -B /work:/work \
    "$CONTAINER" \
    busco \
        -i "$FASTA" \
        -l "$LINEAGE" \
        -m transcriptome \
        -o busco_representatives \
        -c "$CPUS"
