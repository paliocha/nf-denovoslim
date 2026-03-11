#!/bin/bash
#SBATCH --job-name=miniprot_BMED
#SBATCH --output=runs/miniprot_BMED_%j.out
#SBATCH --error=runs/miniprot_BMED_%j.err
#SBATCH --partition=orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00

set -euo pipefail

MINIPROT="$HOME/micromamba/envs/miniprot/bin/miniprot"
OUTDIR="$PROJECTS/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/miniprot_validation"
mkdir -p "$OUTDIR"

GENOME="/mnt/project/FjellheimLab/EGAPx/Briza_maxima/complete.genomic.fna"
PROTEINS="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/BMED/proteins/BMED.faa"
SPECIES="BMED"

echo "=== $SPECIES proteins vs Briza maxima genome ==="
echo "Started: $(date)"
$MINIPROT -t ${SLURM_CPUS_PER_TASK} --gff "$GENOME" "$PROTEINS" > "$OUTDIR/${SPECIES}_vs_BMAX.gff" 2> "$OUTDIR/${SPECIES}_vs_BMAX.log"
echo "Finished miniprot: $(date)"

echo ""
echo "=== Analyzing hit rates by protein length ==="

GFF="$OUTDIR/${SPECIES}_vs_BMAX.gff"

# Get mapped protein IDs (from mRNA lines, Target= field)
awk -F'\t' '$3=="mRNA"{match($9, /Target=([^ ;]+)/, m); if(m[1]) print m[1]}' "$GFF" | sort -u > "$OUTDIR/${SPECIES}_mapped_ids.txt"

# Get all protein IDs and lengths
awk '/^>/{if(id && seq) print id"\t"length(seq); id=substr($1,2); seq=""; next}{seq=seq$0} END{if(id && seq) print id"\t"length(seq)}' "$PROTEINS" > "$OUTDIR/${SPECIES}_protein_lengths.tsv"

TOTAL=$(wc -l < "$OUTDIR/${SPECIES}_protein_lengths.tsv")
MAPPED=$(wc -l < "$OUTDIR/${SPECIES}_mapped_ids.txt")
echo "Total proteins: $TOTAL"
echo "Mapped (any hit): $MAPPED ($(awk "BEGIN{printf \"%.1f\", $MAPPED/$TOTAL*100}")%)"
echo ""

# Join: mark each protein as mapped or not, bin by length
awk -F'\t' '
NR==FNR {mapped[$1]=1; next}
{
    id=$1; len=$2
    hit = (id in mapped) ? 1 : 0

    if      (len < 50)  bin="<50"
    else if (len < 100) bin="50-99"
    else if (len < 150) bin="100-149"
    else if (len < 200) bin="150-199"
    else if (len < 300) bin="200-299"
    else if (len < 500) bin="300-499"
    else                bin="500+"

    total[bin]++
    hits[bin] += hit
    order[bin] = (bin=="<50") ? 1 : (bin=="50-99") ? 2 : (bin=="100-149") ? 3 : (bin=="150-199") ? 4 : (bin=="200-299") ? 5 : (bin=="300-499") ? 6 : 7
}
END {
    printf "%-12s %8s %8s %8s\n", "Length_bin", "Total", "Mapped", "Hit_%"
    printf "%-12s %8s %8s %8s\n", "----------", "-----", "------", "-----"
    for(b in total) rows[order[b]] = sprintf("%-12s %8d %8d %7.1f%%", b, total[b], hits[b], hits[b]/total[b]*100)
    for(i=1; i<=7; i++) if(rows[i]) print rows[i]
}' "$OUTDIR/${SPECIES}_mapped_ids.txt" "$OUTDIR/${SPECIES}_protein_lengths.tsv"

echo ""
echo "=== Done: $(date) ==="
