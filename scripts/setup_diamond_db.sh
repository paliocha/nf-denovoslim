#!/bin/bash
#
# Setup Diamond database for frameshift correction
# Extracts SwissProt FASTA from existing MMseqs2 DB and builds Diamond index
#

set -euo pipefail

# Paths
MMSEQS_DB="/mnt/project/glowberry/transannot/db/SwissProtDB"
OUTPUT_DIR="/mnt/project/glowberry/transannot/db"
DIAMOND_DB="${OUTPUT_DIR}/SwissProtDB_diamond"
FASTA_TMP="${OUTPUT_DIR}/swissprot_tmp.fasta"

# Load modules (if needed on your HPC)
# module load mmseqs2/15.6f452
# module load diamond/2.1.22

echo "=== Diamond DB Setup for Frameshift Correction ==="
echo "Source MMseqs2 DB: ${MMSEQS_DB}"
echo "Output Diamond DB: ${DIAMOND_DB}.dmnd"
echo ""

# Step 1: Extract FASTA from MMseqs2 database
if [ ! -f "${FASTA_TMP}" ]; then
    echo "[1/3] Extracting FASTA from MMseqs2 SwissProt DB..."
    apptainer exec \
        docker://quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_3 \
        mmseqs convert2fasta "${MMSEQS_DB}" "${FASTA_TMP}"
    echo "      → FASTA written to ${FASTA_TMP}"
else
    echo "[1/3] FASTA already exists: ${FASTA_TMP}"
fi

# Step 2: Build Diamond database
if [ ! -f "${DIAMOND_DB}.dmnd" ]; then
    echo "[2/3] Building Diamond database (this takes ~2-5 minutes)..."
    apptainer exec \
        docker://quay.io/biocontainers/diamond:2.1.22--h13889ed_0 \
        diamond makedb \
            --in "${FASTA_TMP}" \
            --db "${DIAMOND_DB}" \
            --threads 8
    echo "      → Diamond DB created: ${DIAMOND_DB}.dmnd"
else
    echo "[2/3] Diamond DB already exists: ${DIAMOND_DB}.dmnd"
fi

# Step 3: Clean up temporary FASTA (optional)
echo "[3/3] Cleaning up temporary FASTA..."
rm -f "${FASTA_TMP}"
echo "      → Cleanup complete"

# Summary
echo ""
echo "=== Setup Complete ==="
echo "Diamond database: ${DIAMOND_DB}.dmnd"
echo ""
echo "Add this to your run scripts (run_BMAX.sh, run_BMED.sh, run_FPRA.sh):"
echo "  --diamond_db ${DIAMOND_DB}.dmnd \\"
echo ""
