#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --job-name=denovoslim-FPRA
#SBATCH --output=%x_%A.out
#SBATCH --chdir=/net/fs-2/scale/OrionStore/Home/martpali/AnnualPerennial/nf-denovoslim/runs/FPRA

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

# Run under fjellheimlab group so all files (work dirs, staging, outputs)
# count against the 40T project quota, not the 500G personal quota.
sg fjellheimlab -c "
nextflow run $PIPELINE_DIR/main.nf \
    -profile apptainer,orion \
    -resume \
    -ansi-log false \
    -w $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/FPRA/work \
    --trinity_fasta $PROJECTS/FjellheimLab/martpali/AnnualPerennial/assemblies/FPRA-Trinity1/FPRA-Trinity.fasta \
    --samplesheet $HOME/AnnualPerennial/nf-denovoslim/FPRA.samplesheet.csv \
    --species_label FPRA \
    --mmseqs2_swissprot $PROJECTS/glowberry/transannot/db/SwissProtDB \
    --mmseqs2_pfam $PROJECTS/glowberry/transannot/db/PfamDB \
    --mmseqs2_eggnog $PROJECTS/glowberry/transannot/db/eggNOG7_profiles \
    --mmseqs2_taxonomy_db $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db/UniRef50taxdb \
    --diamond_db $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db/uniref50.dmnd \
    --eggnog_annotations $PROJECTS/glowberry/transannot/db/e7_as_e5_annotations.tsv \
    --busco_lineage poales_odb12 \
    --filter_taxon 33090 \
    --unix_group fjellheimlab \
    --orion_exclude_nodes cn-37 \
    --outdir $PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/FPRA
"
