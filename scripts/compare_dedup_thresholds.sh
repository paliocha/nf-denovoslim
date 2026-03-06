#!/bin/bash
set -euo pipefail
#
# compare_dedup_thresholds.sh — Post-hoc BUSCO-D comparison at 90% vs 95% nt identity
#
# Runs MMseqs2 clustering at 95% on the same taxonomy-filtered FASTA,
# then BUSCO on both representative sets, and compares Duplicated counts.
#
# Usage:
#   compare_dedup_thresholds.sh \
#       --fasta   <taxonomy-filtered FASTA>  \
#       --species <label, e.g. BMED>         \
#       --lineage <BUSCO lineage, e.g. poales_odb12> \
#       --outdir  <output directory>         \
#       [--cpus 16]                          \
#       [--busco_download_path <path>]
#
# Requires: mmseqs2, busco in PATH (or loaded via module / container).
# Designed to be run standalone after a pipeline run completes.

CPUS=16
BUSCO_DL_PATH=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta)   FASTA="$2";   shift 2 ;;
        --species) SPECIES="$2"; shift 2 ;;
        --lineage) LINEAGE="$2"; shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --cpus)    CPUS="$2";    shift 2 ;;
        --busco_download_path) BUSCO_DL_PATH="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

for var in FASTA SPECIES LINEAGE OUTDIR; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: --$(echo $var | tr 'A-Z' 'a-z') is required" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"

echo "=== BUSCO-D comparison: 90% vs 95% nucleotide identity ==="
echo "Input FASTA:  $FASTA"
echo "Species:      $SPECIES"
echo "BUSCO lineage: $LINEAGE"
echo "Output dir:   $OUTDIR"
echo ""

N_INPUT=$(grep -c '^>' "$FASTA")
echo "Input sequences: $N_INPUT"

# ── Step 1: Cluster at 95% nt identity ──────────────────────────────────

echo ""
echo "[1/4] Clustering at 95% nucleotide identity..."

MMSEQS_DIR="$OUTDIR/mmseqs2_95"
mkdir -p "$MMSEQS_DIR"

mmseqs createdb "$FASTA" "$MMSEQS_DIR/queryDB"
mmseqs cluster "$MMSEQS_DIR/queryDB" "$MMSEQS_DIR/clusterDB" "$MMSEQS_DIR/tmp" \
    --min-seq-id 0.95 \
    -c 0.8 \
    --cov-mode 1 \
    --cluster-mode 2 \
    --threads "$CPUS"
mmseqs createsubdb "$MMSEQS_DIR/clusterDB" "$MMSEQS_DIR/queryDB" "$MMSEQS_DIR/repDB"
mmseqs convert2fasta "$MMSEQS_DIR/repDB" "$OUTDIR/${SPECIES}_reps_95pct.fasta"

N_95=$(grep -c '^>' "$OUTDIR/${SPECIES}_reps_95pct.fasta")
echo "  Representatives at 95%: $N_95"

# ── Step 2: Extract 90% representatives (already the pipeline output) ───

# The input FASTA is the pre-dedup set. We need to also cluster at 90%.
echo ""
echo "[2/4] Clustering at 90% nucleotide identity..."

MMSEQS_DIR_90="$OUTDIR/mmseqs2_90"
mkdir -p "$MMSEQS_DIR_90"

mmseqs createdb "$FASTA" "$MMSEQS_DIR_90/queryDB"
mmseqs cluster "$MMSEQS_DIR_90/queryDB" "$MMSEQS_DIR_90/clusterDB" "$MMSEQS_DIR_90/tmp" \
    --min-seq-id 0.90 \
    -c 0.8 \
    --cov-mode 1 \
    --cluster-mode 2 \
    --threads "$CPUS"
mmseqs createsubdb "$MMSEQS_DIR_90/clusterDB" "$MMSEQS_DIR_90/queryDB" "$MMSEQS_DIR_90/repDB"
mmseqs convert2fasta "$MMSEQS_DIR_90/repDB" "$OUTDIR/${SPECIES}_reps_90pct.fasta"

N_90=$(grep -c '^>' "$OUTDIR/${SPECIES}_reps_90pct.fasta")
echo "  Representatives at 90%: $N_90"

# ── Step 3: BUSCO on both ──────────────────────────────────────────────

BUSCO_DL_ARG=""
if [[ -n "$BUSCO_DL_PATH" ]]; then
    BUSCO_DL_ARG="--download_path $BUSCO_DL_PATH"
fi

echo ""
echo "[3/4] Running BUSCO on 90% representatives..."

busco -i "$OUTDIR/${SPECIES}_reps_90pct.fasta" \
    -o "${SPECIES}_busco_90pct" \
    --out_path "$OUTDIR" \
    -l "$LINEAGE" \
    -m transcriptome \
    -c "$CPUS" \
    $BUSCO_DL_ARG \
    --offline 2>/dev/null || true

echo ""
echo "[4/4] Running BUSCO on 95% representatives..."

busco -i "$OUTDIR/${SPECIES}_reps_95pct.fasta" \
    -o "${SPECIES}_busco_95pct" \
    --out_path "$OUTDIR" \
    -l "$LINEAGE" \
    -m transcriptome \
    -c "$CPUS" \
    $BUSCO_DL_ARG \
    --offline 2>/dev/null || true

# ── Step 4: Compare results ───────────────────────────────────────────

echo ""
echo "=== COMPARISON RESULTS ==="
echo ""

REPORT="$OUTDIR/${SPECIES}_dedup_comparison.txt"

{
    echo "================================================================"
    echo "  BUSCO-D Comparison: 90% vs 95% Nucleotide Identity Dedup"
    echo "  Species: $SPECIES"
    echo "================================================================"
    echo ""
    echo "Input sequences (post-taxonomy): $N_INPUT"
    echo ""
    printf "%-30s %12s %12s\n" "Metric" "90% nt" "95% nt"
    printf "%-30s %12s %12s\n" "------------------------------" "------------" "------------"
    printf "%-30s %12s %12s\n" "Representatives" "$N_90" "$N_95"
    printf "%-30s %12s %12s\n" "Reduction from input" \
        "$(echo "scale=1; ($N_INPUT - $N_90) * 100 / $N_INPUT" | bc)%" \
        "$(echo "scale=1; ($N_INPUT - $N_95) * 100 / $N_INPUT" | bc)%"
    echo ""

    # Parse BUSCO short summaries
    for pct in 90 95; do
        summary=$(find "$OUTDIR/${SPECIES}_busco_${pct}pct" -name "short_summary*.txt" 2>/dev/null | head -1)
        if [[ -n "$summary" && -f "$summary" ]]; then
            echo "--- BUSCO short summary (${pct}% nt identity) ---"
            cat "$summary"
            echo ""
        else
            echo "--- BUSCO summary for ${pct}% not found ---"
            echo ""
        fi
    done

    # Extract key BUSCO numbers for side-by-side comparison
    echo "--- Side-by-side BUSCO comparison ---"
    echo ""
    printf "%-20s %12s %12s\n" "Category" "90% nt" "95% nt"
    printf "%-20s %12s %12s\n" "--------------------" "------------" "------------"

    for pct in 90 95; do
        summary=$(find "$OUTDIR/${SPECIES}_busco_${pct}pct" -name "short_summary*.txt" 2>/dev/null | head -1)
        if [[ -n "$summary" && -f "$summary" ]]; then
            eval "C_${pct}=$(grep -oP '\d+(?=\s+Complete BUSCOs)' "$summary" 2>/dev/null || echo 'N/A')"
            eval "S_${pct}=$(grep -oP '\d+(?=\s+Complete and single-copy)' "$summary" 2>/dev/null || echo 'N/A')"
            eval "D_${pct}=$(grep -oP '\d+(?=\s+Complete and duplicated)' "$summary" 2>/dev/null || echo 'N/A')"
            eval "F_${pct}=$(grep -oP '\d+(?=\s+Fragmented)' "$summary" 2>/dev/null || echo 'N/A')"
            eval "M_${pct}=$(grep -oP '\d+(?=\s+Missing)' "$summary" 2>/dev/null || echo 'N/A')"
        fi
    done

    printf "%-20s %12s %12s\n" "Complete (C)" "${C_90:-N/A}" "${C_95:-N/A}"
    printf "%-20s %12s %12s\n" "Single-copy (S)" "${S_90:-N/A}" "${S_95:-N/A}"
    printf "%-20s %12s %12s\n" "Duplicated (D)" "${D_90:-N/A}" "${D_95:-N/A}"
    printf "%-20s %12s %12s\n" "Fragmented (F)" "${F_90:-N/A}" "${F_95:-N/A}"
    printf "%-20s %12s %12s\n" "Missing (M)" "${M_90:-N/A}" "${M_95:-N/A}"
    echo ""
    echo "Key metric: D count decrease from 95% → 90% shows how many"
    echo "near-duplicate BUSCOs the stricter threshold collapses."
    echo "If M (Missing) increases, the threshold is too aggressive."

} | tee "$REPORT"

echo ""
echo "Full report written to: $REPORT"
