#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=download-pfam
#SBATCH --output=%x_%A.out

# Download Pfam-A.hmm for the HMMER_EXTEND protein extension step.
#
# Pfam-A.hmm contains ~20K protein family HMM profiles (~1.5 GB compressed).
# pyhmmer reads the plain .hmm file directly — no hmmpress needed.

set -euo pipefail

DB_DIR="${PROJECTS}/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

echo "Downloading Pfam-A.hmm.gz from InterPro/EBI..."
wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

echo "Decompressing..."
gunzip -f Pfam-A.hmm.gz

echo "Done. Files:"
ls -lh Pfam-A.hmm*

echo ""
echo "Pfam-A.hmm is ready at: ${DB_DIR}/Pfam-A.hmm"
echo "Use --pfam_hmm ${DB_DIR}/Pfam-A.hmm in your pipeline run."
