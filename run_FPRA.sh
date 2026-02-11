#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --job-name=denovoslim-FPRA
#SBATCH --output=denovoslim-FPRA_%A.out

source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate Nextflow

module load Java
module load Anaconda3
module load singularity

PIPELINE_DIR=$HOME/AnnualPerennial/nf-denovoslim
LAUNCH_DIR=$HOME/AnnualPerennial/nf-denovoslim/runs/FPRA
mkdir -p $LAUNCH_DIR
cd $LAUNCH_DIR

nextflow run $PIPELINE_DIR/main.nf \
    -profile apptainer,slurm \
    -resume \
    -ansi-log false \
    -w $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/FPRA/work \
    --trinity_fasta $PROJECTS/FjellheimLab/martpali/AnnualPerennial/assemblies/FPRA-Trinity1/FPRA-Trinity.fasta \
    --samplesheet $HOME/AnnualPerennial/nf-denovoslim/FPRA.samplesheet.csv \
    --species_label FPRA \
    --mmseqs2_swissprot $PROJECTS/glowberry/transannot/db/SwissProtDB \
    --mmseqs2_pfam $PROJECTS/glowberry/transannot/db/PfamDB \
    --mmseqs2_eggnog $PROJECTS/glowberry/transannot/db/eggNOG7_profiles \
    --outdir $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/FPRA
