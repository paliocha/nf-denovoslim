#!/bin/bash
#SBATCH --partition=orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --job-name=download_uniref50
#SBATCH --output=download_uniref50_%j.out
#SBATCH --error=download_uniref50_%j.err

##
## Download all UniRef50 prerequisites and then build DIAMOND + MMseqs2
## taxonomy databases by submitting the two build scripts as SLURM jobs.
##
## Usage:
##   sbatch scripts/download_uniref50.sh
##
## This script downloads:
##   1. UniRef50 FASTA           (~12 GB compressed)
##   2. NCBI taxonomy dump       (~70 MB compressed)
##   3. UniProt ID-mapping file  (~100 GB compressed)
##
## After downloads complete, it derives the taxidmapping_prefixed file
## and then submits:
##   - scripts/build_uniref50_diamond.sh      (DIAMOND makedb)
##   - scripts/build_uniref50_mmseqstaxdb.sh  (MMseqs2 createdb + createtaxdb)
##
## On Orion, compute and login nodes CAN reach the internet.
##

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

DB_DIR="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db"
TMP_DIR="$DB_DIR/tmp_uniref50"
mkdir -p "$DB_DIR" "$TMP_DIR"

echo "=================================================="
echo "UniRef50 database download + build orchestrator"
echo "Started: $(date)"
echo "DB dir:  $DB_DIR"
echo "Tmp dir: $TMP_DIR"
echo "=================================================="
echo ""

# ── 1. Download UniRef50 FASTA ──────────────────────────────────────────
cd "$DB_DIR"
if [ -f "uniref50.fasta.gz" ]; then
    echo "[$(date)] uniref50.fasta.gz already exists ($(ls -lh uniref50.fasta.gz | awk '{print $5}')) — skipping download"
else
    echo "[$(date)] Downloading UniRef50 FASTA (~12 GB)..."
    wget -c -q --show-progress \
        https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
    echo "[$(date)] UniRef50 FASTA download complete: $(ls -lh uniref50.fasta.gz | awk '{print $5}')"
fi
echo ""

# ── 2. Download NCBI taxonomy dump ──────────────────────────────────────
cd "$TMP_DIR"
if [ -f "nodes.dmp" ] && [ -f "names.dmp" ] && [ -f "merged.dmp" ] && [ -f "delnodes.dmp" ]; then
    echo "[$(date)] NCBI taxonomy dump already present — skipping download"
else
    echo "[$(date)] Downloading NCBI taxonomy dump..."
    wget -O taxdump.tar.gz \
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf taxdump.tar.gz names.dmp nodes.dmp merged.dmp delnodes.dmp
    echo "[$(date)] Taxonomy dump extracted: names.dmp nodes.dmp merged.dmp delnodes.dmp"
fi
echo ""

# ── 3. Download UniProt ID-mapping and derive tax mapping ──────────────
if [ -f "taxidmapping_prefixed" ]; then
    echo "[$(date)] taxidmapping_prefixed already exists ($(wc -l < taxidmapping_prefixed) lines) — skipping"
else
    if [ ! -f "idmapping.dat.gz" ]; then
        echo "[$(date)] Downloading UniProt idmapping.dat.gz (~100 GB)..."
        echo "          This is the largest download and will take a while."
        wget -c -q --show-progress -O idmapping.dat.gz \
            https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
        echo "[$(date)] idmapping.dat.gz download complete: $(ls -lh idmapping.dat.gz | awk '{print $5}')"
    else
        echo "[$(date)] idmapping.dat.gz already exists — skipping download"
    fi

    echo "[$(date)] Extracting NCBI_TaxID entries with UniRef50_ prefix..."
    gunzip -c idmapping.dat.gz \
        | awk '$2 == "NCBI_TaxID" {print "UniRef50_" $1 "\t" $3}' \
        > taxidmapping_prefixed
    echo "[$(date)] taxidmapping_prefixed created: $(wc -l < taxidmapping_prefixed) entries"
fi
echo ""

# ── 4. Validate all prerequisites ──────────────────────────────────────
echo "=== Prerequisite check ==="
MISSING=0
for f in "$DB_DIR/uniref50.fasta.gz" \
         "$TMP_DIR/names.dmp" "$TMP_DIR/nodes.dmp" \
         "$TMP_DIR/merged.dmp" "$TMP_DIR/delnodes.dmp" \
         "$TMP_DIR/taxidmapping_prefixed"; do
    if [ -f "$f" ]; then
        echo "  OK: $(basename $f)"
    else
        echo "  MISSING: $f"
        MISSING=1
    fi
done

if [ "$MISSING" -eq 1 ]; then
    echo ""
    echo "ERROR: Some prerequisite files are missing. Cannot submit build jobs."
    exit 1
fi
echo ""

# ── 5. Submit build jobs ───────────────────────────────────────────────
echo "=== Submitting build jobs ==="

DIAMOND_JOB=$(sbatch --parsable "$PIPELINE_DIR/scripts/build_uniref50_diamond.sh")
echo "  DIAMOND build submitted: job $DIAMOND_JOB"

MMSEQS_JOB=$(sbatch --parsable "$PIPELINE_DIR/scripts/build_uniref50_mmseqstaxdb.sh")
echo "  MMseqs2 taxdb build submitted: job $MMSEQS_JOB"

echo ""
echo "=== All downloads complete, build jobs submitted ==="
echo "Finished: $(date)"
echo ""
echo "Monitor builds with:"
echo "  sacct -j $DIAMOND_JOB,$MMSEQS_JOB --format=JobID,JobName%30,State,Elapsed"
echo ""
echo "Once both complete, the databases are ready:"
echo "  --diamond_db          $DB_DIR/uniref50.dmnd"
echo "  --mmseqs2_taxonomy_db $DB_DIR/UniRef50taxdb"
