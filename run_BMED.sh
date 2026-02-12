#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --job-name=denovoslim-BMED
#SBATCH --output=denovoslim-BMED_%A.out

source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate Nextflow

module load Java
module load Anaconda3
module load singularity

PIPELINE_DIR=$HOME/AnnualPerennial/nf-denovoslim
LAUNCH_DIR=$HOME/AnnualPerennial/nf-denovoslim/runs/BMED
mkdir -p $LAUNCH_DIR
cd $LAUNCH_DIR

nextflow run $PIPELINE_DIR/main.nf \
    -profile apptainer,slurm \
    -resume \
    -ansi-log false \
    -w $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/BMED/work \
    --trinity_fasta $PROJECTS/FjellheimLab/martpali/AnnualPerennial/assemblies/BMED-Trinity1/BMED-Trinity.fasta \
    --samplesheet $HOME/AnnualPerennial/nf-denovoslim/BMED.samplesheet.csv \
    --species_label BMED \
    --mmseqs2_swissprot $PROJECTS/glowberry/transannot/db/SwissProtDB \
    --mmseqs2_pfam $PROJECTS/glowberry/transannot/db/PfamDB \
    --mmseqs2_eggnog $PROJECTS/glowberry/transannot/db/eggNOG7_profiles \
    --mmseqs2_taxonomy_db $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db/UniProtTrEMBLtaxdb \
    --diamond_db $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db/uniref90.dmnd \
    --filter_taxon 35493 \
    --outdir $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/BMED
