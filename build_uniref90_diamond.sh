#!/bin/bash
#SBATCH --job-name=diamond_uniref90
#SBATCH --partition=orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=diamond_uniref90_%j.out
#SBATCH --error=diamond_uniref90_%j.err

# Build UniRef90 Diamond database for frameshift correction
# Expected: ~30-40 GB compressed FASTA download, ~50-80 GB final .dmnd

set -euo pipefail

DB_DIR="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

echo "=================================================="
echo "Building UniRef90 Diamond database"
echo "Started: $(date)"
echo "Output: ${DB_DIR}/uniref90.dmnd"
echo "=================================================="
echo ""

# 1. Download UniRef90 from UniProt
if [ ! -f uniref90.fasta.gz ]; then
    echo "[$(date)] Downloading UniRef90 FASTA (~30-40 GB)..."
    wget -c https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
    echo "[$(date)] Download complete."
else
    echo "[$(date)] UniRef90 FASTA already downloaded, skipping."
fi

# 2. Build Diamond database
echo ""
echo "[$(date)] Building Diamond database (this will take 30-60 minutes)..."
apptainer exec docker://quay.io/biocontainers/diamond:2.1.22--h13889ed_0 \
    diamond makedb \
        --in uniref90.fasta.gz \
        -d uniref90 \
        -p ${SLURM_CPUS_PER_TASK:-16}

echo ""
echo "[$(date)] Diamond database build complete!"
echo "Database: ${DB_DIR}/uniref90.dmnd"
echo "Size: $(du -h uniref90.dmnd | cut -f1)"
echo ""
echo "Add to your run scripts:"
echo "  --diamond_db ${DB_DIR}/uniref90.dmnd"
echo ""
echo "=================================================="
echo "Finished: $(date)"
echo "=================================================="
