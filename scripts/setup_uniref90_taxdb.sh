#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --job-name=build-uniref90-taxdb
#SBATCH --output=%x_%A.out

##
## Build a UniRef90 MMseqs2 taxonomy database (seqTaxDB) from scratch.
##
## UniRef90 has ~90M clusters (vs ~252M in full TrEMBL), giving a ~3.5× smaller
## k-mer index that fits comfortably in 300–400 GB RAM without target splitting.
## This dramatically speeds up the MMseqs2 taxonomy step.
##
## Usage:
##   sbatch scripts/setup_uniref90_taxdb.sh
##
## The database will be created at:
##   $DB_DIR/UniRef90taxdb
##

set -euo pipefail

# ── Configuration ──
DB_DIR="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db"
DB_NAME="UniRef90taxdb"
THREADS=32

# Load modules / containers
module load singularity 2>/dev/null || true
MMSEQS_IMG="quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0"

# Use apptainer exec if available, otherwise assume mmseqs is in PATH
if command -v apptainer &>/dev/null; then
    MMSEQS="apptainer exec --bind /mnt/project:/mnt/project --bind /work:/work $MMSEQS_IMG mmseqs"
elif command -v singularity &>/dev/null; then
    MMSEQS="singularity exec --bind /mnt/project:/mnt/project --bind /work:/work $MMSEQS_IMG mmseqs"
else
    MMSEQS="mmseqs"
fi

cd "$DB_DIR"
mkdir -p tmp_uniref90

echo "=== Downloading UniRef90 via mmseqs databases ==="
echo "This downloads the FASTA, taxonomy mapping, and NCBI taxonomy dump."
echo "Started: $(date)"

$MMSEQS databases UniRef90 "$DB_NAME" tmp_uniref90 --threads $THREADS

echo ""
echo "=== Download & DB creation complete ==="
echo "Finished: $(date)"
echo ""
echo "Database files:"
ls -lh ${DB_NAME}*
echo ""
echo "Sequence count:"
grep -c '.' ${DB_NAME}.index
echo ""
echo "=== Done ==="
echo "Update your run scripts to use:"
echo "  --mmseqs2_taxonomy_db $DB_DIR/$DB_NAME"

# Cleanup
rm -rf tmp_uniref90
