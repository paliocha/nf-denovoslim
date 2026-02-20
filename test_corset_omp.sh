#!/bin/bash
#SBATCH --job-name=corset_FPRA
#SBATCH --partition=orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=corset_FPRA_%j.out
#SBATCH --error=corset_FPRA_%j.err
set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────
SIF="/mnt/users/martpali/AnnualPerennial/nf-denovoslim/containers/corset/corset_1.10.sif"
WORK="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/FPRA/work"
OUTDIR="$TMPDIR/corset_test"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ── Stage eq_classes via symlinked directory trees ───────────────────
# Corset expects  <sample>_quant/aux_info/eq_classes.txt
SAMPLES=(
  FPRA03_T1_L FPRA05_T1_L FPRA07_T1_L FPRA08_T1_L
  FPRA03_T1_R FPRA05_T1_R FPRA07_T1_R FPRA08_T1_R
  FPRA16_T2_L FPRA18_T2_L FPRA19_T2_L FPRA21_T2_L
  FPRA16_T2_R FPRA18_T2_R FPRA19_T2_R FPRA21_T2_R
  FPRA35_T3_L FPRA36_T3_L FPRA37_T3_L FPRA38_T3_L
  FPRA35_T3_R FPRA36_T3_R FPRA37_T3_R FPRA38_T3_R
  FPRA55_T4_L FPRA56_T4_L FPRA58_T4_L FPRA60_T4_L
  FPRA55_T4_R FPRA56_T4_R FPRA58_T4_R FPRA60_T4_R
  FPRA75_T5_L FPRA76_T5_L FPRA78_T5_L FPRA80_T5_L
  FPRA75_T5_R FPRA76_T5_R FPRA78_T5_R FPRA80_T5_R
)

# Real paths to eq_classes.txt for each sample (same order as SAMPLES)
EQPATHS=(
  "$WORK/80/df07120d901beb78af467fe1ef9696/FPRA03_T1_L_quant/aux_info/eq_classes.txt"
  "$WORK/71/19ca8076d232ab8e4e4529d1f55a54/FPRA05_T1_L_quant/aux_info/eq_classes.txt"
  "$WORK/ac/2d52560a87bac9e439b0e0be316fc1/FPRA07_T1_L_quant/aux_info/eq_classes.txt"
  "$WORK/ff/6fb94eb5ba6aa7e5744234af3fda9f/FPRA08_T1_L_quant/aux_info/eq_classes.txt"
  "$WORK/3c/9b5be17eb20853e8245dd876de115c/FPRA03_T1_R_quant/aux_info/eq_classes.txt"
  "$WORK/9c/7058761b849577e9a7619ee460a4d3/FPRA05_T1_R_quant/aux_info/eq_classes.txt"
  "$WORK/e5/4e34ea4a61a0d6a9a9ef29c0b426a3/FPRA07_T1_R_quant/aux_info/eq_classes.txt"
  "$WORK/90/e11865d45a753bffd4350aae3f7093/FPRA08_T1_R_quant/aux_info/eq_classes.txt"
  "$WORK/29/0bac2b1d715ab03e7239e020216550/FPRA16_T2_L_quant/aux_info/eq_classes.txt"
  "$WORK/e0/bfadc8e4536896ed6fa0901c4a4867/FPRA18_T2_L_quant/aux_info/eq_classes.txt"
  "$WORK/bb/d2bf17752679d85b4f8f09ecf7ea2e/FPRA19_T2_L_quant/aux_info/eq_classes.txt"
  "$WORK/56/fe8f2a30cf9ad95f695c7ad73d5223/FPRA21_T2_L_quant/aux_info/eq_classes.txt"
  "$WORK/b4/023abd96c71ae35947157852fe1045/FPRA16_T2_R_quant/aux_info/eq_classes.txt"
  "$WORK/3e/4762a5323c1a66a8c56c6deb955215/FPRA18_T2_R_quant/aux_info/eq_classes.txt"
  "$WORK/1c/869ffe41dab489f5ed14cb0de9be68/FPRA19_T2_R_quant/aux_info/eq_classes.txt"
  "$WORK/0b/a94030eb21c842f5d8f88271a8bd91/FPRA21_T2_R_quant/aux_info/eq_classes.txt"
  "$WORK/d8/da25025348c2ca955d33e23c1356ff/FPRA35_T3_L_quant/aux_info/eq_classes.txt"
  "$WORK/60/9ba1ba125919e44d18a51cb7f2a8eb/FPRA36_T3_L_quant/aux_info/eq_classes.txt"
  "$WORK/42/db2011d0711d25c06d6d088bfecaeb/FPRA37_T3_L_quant/aux_info/eq_classes.txt"
  "$WORK/b4/8b66523102157034ee8e858081ef90/FPRA38_T3_L_quant/aux_info/eq_classes.txt"
  "$WORK/6b/2053e47c8a9496673aadbc636cad3e/FPRA35_T3_R_quant/aux_info/eq_classes.txt"
  "$WORK/7b/45bedf43d446de32ee7c194d93b3a3/FPRA36_T3_R_quant/aux_info/eq_classes.txt"
  "$WORK/4d/7a07fcb5ffc31ded05b286a6f6aab1/FPRA37_T3_R_quant/aux_info/eq_classes.txt"
  "$WORK/21/317298864d50ea851964399a0d26ec/FPRA38_T3_R_quant/aux_info/eq_classes.txt"
  "$WORK/74/c9d89205ef7fcb675a3f267d1b82ab/FPRA55_T4_L_quant/aux_info/eq_classes.txt"
  "$WORK/76/06d02eb322af58609c450ac9176e30/FPRA56_T4_L_quant/aux_info/eq_classes.txt"
  "$WORK/28/2dbecb82c14633a7c249ce94facd77/FPRA58_T4_L_quant/aux_info/eq_classes.txt"
  "$WORK/21/63df4b2eaaa0562acc6555e08940fd/FPRA60_T4_L_quant/aux_info/eq_classes.txt"
  "$WORK/04/734d0564b3aeeb38531c6c0ed5d9ae/FPRA55_T4_R_quant/aux_info/eq_classes.txt"
  "$WORK/c4/3ba373065d024d62ea577b40a7c6f9/FPRA56_T4_R_quant/aux_info/eq_classes.txt"
  "$WORK/ad/cf8b93c6656402f55330cffccc11b8/FPRA58_T4_R_quant/aux_info/eq_classes.txt"
  "$WORK/75/f4213ba02880d7b52063bd5e03707d/FPRA60_T4_R_quant/aux_info/eq_classes.txt"
  "$WORK/83/4c05a4caf60cd90086d792084729fc/FPRA75_T5_L_quant/aux_info/eq_classes.txt"
  "$WORK/f4/6982a656e6d4209c1cb455a8c868ba/FPRA76_T5_L_quant/aux_info/eq_classes.txt"
  "$WORK/57/d515d2cc4ed2f17759a046c151b5a8/FPRA78_T5_L_quant/aux_info/eq_classes.txt"
  "$WORK/11/b362e5569cb60cc8f9ebd565d94dd1/FPRA80_T5_L_quant/aux_info/eq_classes.txt"
  "$WORK/e4/f0249c84317b5122e5074b56411db3/FPRA75_T5_R_quant/aux_info/eq_classes.txt"
  "$WORK/0c/0a6b365bfcd94b4a4d06e370fd17c9/FPRA76_T5_R_quant/aux_info/eq_classes.txt"
  "$WORK/69/2a98815de33ccc92a8c5a7593355d3/FPRA78_T5_R_quant/aux_info/eq_classes.txt"
  "$WORK/ba/8263d61239c1117e3d43bb3379aad8/FPRA80_T5_R_quant/aux_info/eq_classes.txt"
)

echo "=== Staging eq_classes symlinks ==="
for i in "${!SAMPLES[@]}"; do
  s="${SAMPLES[$i]}"
  p="${EQPATHS[$i]}"
  mkdir -p "${s}_quant/aux_info"
  ln -sf "$p" "${s}_quant/aux_info/eq_classes.txt"
done
echo "Staged ${#SAMPLES[@]} samples"

# Verify all files are readable
for i in "${!SAMPLES[@]}"; do
  f="${SAMPLES[$i]}_quant/aux_info/eq_classes.txt"
  if [[ ! -r "$f" ]]; then
    echo "ERROR: missing $f" >&2; exit 1
  fi
done
echo "All eq_classes verified"

# ── Build file-argument list (same order) ────────────────────────────
EQ_ARGS=()
for s in "${SAMPLES[@]}"; do
  EQ_ARGS+=("${s}_quant/aux_info/eq_classes.txt")
done

# ── Groups & names (identical to pipeline) ───────────────────────────
# N.B. Cannot use GROUPS — it is a bash built-in (user GID array)
SAMPLE_GROUPS="T1_L,T1_L,T1_L,T1_L,T1_R,T1_R,T1_R,T1_R"
SAMPLE_GROUPS+=",T2_L,T2_L,T2_L,T2_L,T2_R,T2_R,T2_R,T2_R"
SAMPLE_GROUPS+=",T3_L,T3_L,T3_L,T3_L,T3_R,T3_R,T3_R,T3_R"
SAMPLE_GROUPS+=",T4_L,T4_L,T4_L,T4_L,T4_R,T4_R,T4_R,T4_R"
SAMPLE_GROUPS+=",T5_L,T5_L,T5_L,T5_L,T5_R,T5_R,T5_R,T5_R"

NAMES=$(IFS=,; echo "${SAMPLES[*]}")

# ── Run corset-omp ───────────────────────────────────────────────────
echo "=== Running corset with ${SLURM_CPUS_PER_TASK} threads ==="
echo "Start: $(date '+%Y-%m-%d %H:%M:%S')"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

apptainer exec \
    --bind "$OUTDIR:$OUTDIR" \
    --bind "$WORK:$WORK" \
    --env OMP_NUM_THREADS="${OMP_NUM_THREADS}" \
    "$SIF" \
    corset \
        -i salmon_eq_classes \
        -g "$SAMPLE_GROUPS" \
        -n "$NAMES" \
        -p corset_fpra \
        "${EQ_ARGS[@]}"

echo "End:   $(date '+%Y-%m-%d %H:%M:%S')"

# ── Copy results back to persistent storage ──────────────────────────
RESULT_DIR="/mnt/users/martpali/AnnualPerennial/nf-denovoslim/test_corset_results"
mkdir -p "$RESULT_DIR"
cp -v corset_fpra-clusters.txt corset_fpra-counts.txt "$RESULT_DIR/"

# ── Compare with original if available ───────────────────────────────
ORIG_CLUSTERS="$WORK/ae/4e0283314e6241eff5bfed7d57f516/corset-clusters.txt"
if [[ -f "$ORIG_CLUSTERS" ]]; then
  echo ""
  echo "=== Comparison with original corset output ==="
  echo "Original clusters: $(wc -l < "$ORIG_CLUSTERS")"
  echo "New clusters:      $(wc -l < corset_fpra-clusters.txt)"
  echo ""
  # Check if cluster assignments match (column 2)
  diff <(sort "$ORIG_CLUSTERS") <(sort corset_fpra-clusters.txt) > /dev/null 2>&1 \
    && echo "RESULT: Cluster assignments IDENTICAL" \
    || echo "RESULT: Cluster assignments DIFFER (expected — different hash seeding)"
  echo ""
  # Count unique clusters
  echo "Original unique clusters: $(cut -f2 "$ORIG_CLUSTERS" | sort -u | wc -l)"
  echo "New unique clusters:      $(cut -f2 corset_fpra-clusters.txt | sort -u | wc -l)"
fi

echo ""
echo "=== Done ==="
