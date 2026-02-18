#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=2:00:00
#SBATCH --job-name=build-mapping
#SBATCH --output=build-mapping_%A.out

# Build the UniRef90taxdb_mapping file using sort+join.
# The awk-based approach in mmseqs createtaxdb OOMs because it loads
# 137M taxidmapping entries into a hash table.

set -euo pipefail

DB_DIR="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/db"
TMP_DIR="$DB_DIR/tmp_uniref90"
DB_NAME="UniRef90taxdb"

cd "$DB_DIR"

echo "[$(date)] Building _mapping file via sort+join..."

# Extract accession (col 2) and internal ID (col 1) from lookup, sort by accession
echo "[$(date)] Sorting lookup file (188M entries)..."
awk -F'\t' '{print $2"\t"$1}' "${DB_NAME}.lookup" \
    | LC_ALL=C sort -k1,1 -t'	' --parallel=8 -S 80G \
    > "$TMP_DIR/lookup_sorted.tmp"

# Sort the prefixed taxidmapping by accession
echo "[$(date)] Sorting taxidmapping (137M entries)..."
LC_ALL=C sort -k1,1 -t'	' --parallel=8 -S 80G \
    "$TMP_DIR/taxidmapping_prefixed" \
    > "$TMP_DIR/taxidmapping_sorted.tmp"

# Join on accession, output internal_id \t taxid
echo "[$(date)] Joining..."
LC_ALL=C join -t'	' -1 1 -2 1 -o '1.2,2.2' \
    "$TMP_DIR/lookup_sorted.tmp" \
    "$TMP_DIR/taxidmapping_sorted.tmp" \
    > "${DB_NAME}_mapping"

echo "[$(date)] Done. Mapping file:"
wc -l "${DB_NAME}_mapping"
head -5 "${DB_NAME}_mapping"

# Cleanup temp files
rm -f "$TMP_DIR/lookup_sorted.tmp" "$TMP_DIR/taxidmapping_sorted.tmp"

echo "[$(date)] Complete!"
