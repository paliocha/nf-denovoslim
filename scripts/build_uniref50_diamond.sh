#!/bin/bash
#SBATCH --job-name=diamond_uniref50
#SBATCH --partition=orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=diamond_uniref50_%j.out
#SBATCH --error=diamond_uniref50_%j.err

##
## Build UniRef50 Diamond database for frameshift correction.
##
## UniRef50 is ~3-4× smaller than UniRef90 (~52M vs 188M sequences),
## giving proportionally faster DIAMOND blastx searches.  For frameshift
## correction we only need the top hit, so UniRef50 clustering has
## negligible effect on correction accuracy.
##
## Download the FASTA on the login node first:
##
##   cd $DB_DIR
##   wget -c https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
##
## Then:  sbatch scripts/build_uniref50_diamond.sh
##
## Expected output: ~25-30 GB .dmnd file (vs ~86 GB for UniRef90)
##

set -euo pipefail

DB_DIR="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

echo "=================================================="
echo "Building UniRef50 Diamond database"
echo "Started: $(date)"
echo "Output: ${DB_DIR}/uniref50.dmnd"
echo "=================================================="
echo ""

# 1. Check FASTA exists
if [ ! -f uniref50.fasta.gz ]; then
    echo "ERROR: uniref50.fasta.gz not found in $DB_DIR"
    echo "Download it on the login node:"
    echo "  wget -c https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
    exit 1
fi

# 2. Build Diamond database
echo "[$(date)] Building Diamond database (estimated 15-30 minutes)..."
apptainer exec docker://quay.io/biocontainers/diamond:2.1.22--h13889ed_0 \
    diamond makedb \
        --in uniref50.fasta.gz \
        -d uniref50 \
        --threads 16

echo ""
echo "[$(date)] Build complete."
echo ""

# 3. Validate
echo "Database file:"
ls -lh uniref50.dmnd
echo ""
echo "To use, update run scripts:"
echo "  --diamond_db $DB_DIR/uniref50.dmnd"
echo ""
echo "Done."
