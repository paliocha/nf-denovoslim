#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=600G
#SBATCH --time=12:00:00
#SBATCH --job-name=build-uniref90-taxdb
#SBATCH --output=%x_%A.out

##
## Build a UniRef90 MMseqs2 taxonomy database (seqTaxDB) from scratch.
##
## IMPORTANT: Orion compute nodes cannot reach the internet.  Download the
## required files ON THE LOGIN NODE first, then sbatch this script:
##
##   # 1. Download on login node:
##   cd $DB_DIR/tmp_uniref90
##   wget -O taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
##   tar -xzf taxdump.tar.gz names.dmp nodes.dmp merged.dmp delnodes.dmp
##   wget -O idmapping.dat.gz https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
##   gunzip -c idmapping.dat.gz | awk '$2 == "NCBI_TaxID" {print "UniRef90_" $1 "\t" $3}' > taxidmapping_prefixed
##
##   # 2. Submit build job:
##   sbatch scripts/setup_uniref90_taxdb.sh
##
## The database will be created at:
##   $DB_DIR/UniRef90taxdb
##

set -euo pipefail

# ── Configuration ──
DB_DIR="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db"
DB_NAME="UniRef90taxdb"
TMP_DIR="$DB_DIR/tmp_uniref90"
THREADS=32

# Load modules / containers
module load singularity 2>/dev/null || true
MMSEQS_IMG="docker://quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0"

# Use apptainer exec if available, otherwise assume mmseqs is in PATH
if command -v apptainer &>/dev/null; then
    MMSEQS="apptainer exec --bind /mnt/project:/mnt/project --bind /work:/work $MMSEQS_IMG mmseqs"
elif command -v singularity &>/dev/null; then
    MMSEQS="singularity exec --bind /mnt/project:/mnt/project --bind /work:/work $MMSEQS_IMG mmseqs"
else
    MMSEQS="mmseqs"
fi

cd "$DB_DIR"
mkdir -p "$TMP_DIR"

echo "=== Building UniRef90 MMseqs2 taxonomy DB ==="
echo "Started: $(date)"

# Step 1: Create MMseqs2 sequence DB from existing FASTA (skip if already done)
if [ -f "${DB_NAME}.index" ] && [ -f "${DB_NAME}" ]; then
    echo "[$(date)] createdb already complete — skipping ($(wc -l < ${DB_NAME}.index) sequences)"
else
    echo "[$(date)] createdb from uniref90.fasta.gz..."
    $MMSEQS createdb uniref90.fasta.gz "$DB_NAME" --threads "$THREADS"
fi

# Step 2: Verify pre-downloaded files exist (compute nodes have no internet)
echo "[$(date)] Checking pre-downloaded files..."
for f in "$TMP_DIR/names.dmp" "$TMP_DIR/nodes.dmp" "$TMP_DIR/merged.dmp" "$TMP_DIR/taxidmapping_prefixed"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file missing: $f"
        echo "Download files on login node first — see instructions at top of this script."
        exit 1
    fi
done
echo "All pre-downloaded files present."

# Step 3: Remove stale taxonomy outputs so createtaxdb regenerates them
echo "[$(date)] Removing stale taxonomy output files..."
rm -f "${DB_NAME}_mapping" "${DB_NAME}_taxonomy" \
      "${DB_NAME}_names.dmp" "${DB_NAME}_nodes.dmp" "${DB_NAME}_merged.dmp"

# Step 4: Add taxonomy using pre-downloaded files (mmseqs built-in createtaxdb)
echo "[$(date)] createtaxdb (using pre-downloaded taxdump + idmapping)..."
$MMSEQS createtaxdb "$DB_NAME" "$TMP_DIR" --threads "$THREADS" \
    --ncbi-tax-dump "$TMP_DIR" \
    --tax-mapping-file "$TMP_DIR/taxidmapping_prefixed"

# Validate output
echo ""
echo "=== Validation ==="
if [ -s "${DB_NAME}_mapping" ]; then
    echo "_mapping: $(wc -l < "${DB_NAME}_mapping") entries"
    echo "First 5 entries:"
    head -5 "${DB_NAME}_mapping"
else
    echo "WARNING: _mapping file is empty or missing!"
fi

echo ""
echo "=== DB creation complete ==="
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
