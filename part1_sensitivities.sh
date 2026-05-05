#!/usr/bin/env bash
# part1_sensitivities.sh
# Run all five Part 1 sensitivity analyses and regenerate figures.
# Run from the repo root (scotland/) with the PhD conda environment active:
#   conda activate PhD
#   bash part1_sensitivities.sh
# Or without activating first:
#   conda run -n PhD bash part1_sensitivities.sh

set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "======================================================="
echo "  Part 1 sensitivity analyses"
echo "  $(date)"
echo "======================================================="

# ------------------------------------------------------------------
# 1. Health-board spatial clustering
# ------------------------------------------------------------------
echo ""
echo "--- 1/5  --cluster-by health_board ---"
conda run -n PhD python part1/main/main_analysis.py \
    --cluster-by health_board \
    --tables-dir  part1/main/tables_health_board \
    --figures-dir part1/main/figures_health_board \
    --cache-dir   part1/main/cache_health_board

echo "--- figures for health_board ---"
conda run -n PhD python part1/main/manuscript/make_figures.py \
    --tables-dir part1/main/tables_health_board \
    --out-dir    part1/main/manuscript/figures_health_board \
    --cache-dir  part1/main/cache_health_board

# ------------------------------------------------------------------
# 2. Sampling-pool size offset (wn_no_sequences)
# ------------------------------------------------------------------
echo ""
echo "--- 2/5  --use-size-offset ---"
conda run -n PhD python part1/main/main_analysis.py \
    --use-size-offset \
    --tables-dir  part1/main/tables_size_offset \
    --figures-dir part1/main/figures_size_offset \
    --cache-dir   part1/main/cache_size_offset

echo "--- figures for size_offset ---"
conda run -n PhD python part1/main/manuscript/make_figures.py \
    --tables-dir part1/main/tables_size_offset \
    --out-dir    part1/main/manuscript/figures_size_offset \
    --cache-dir  part1/main/cache_size_offset

# ------------------------------------------------------------------
# 3. Index-case SIMD instead of mean-cluster SIMD
# ------------------------------------------------------------------
echo ""
echo "--- 3/5  --use-index-simd ---"
conda run -n PhD python part1/main/main_analysis.py \
    --use-index-simd \
    --tables-dir  part1/main/tables_index_simd \
    --figures-dir part1/main/figures_index_simd \
    --cache-dir   part1/main/cache_index_simd

echo "--- figures for index_simd ---"
conda run -n PhD python part1/main/manuscript/make_figures.py \
    --tables-dir part1/main/tables_index_simd \
    --out-dir    part1/main/manuscript/figures_index_simd \
    --cache-dir  part1/main/cache_index_simd

# ------------------------------------------------------------------
# 4. Tail winsorisation at 99th percentile
# ------------------------------------------------------------------
echo ""
echo "--- 4/5  --winsorise-quantile 0.99 ---"
conda run -n PhD python part1/main/main_analysis.py \
    --winsorise-quantile 0.99 \
    --tables-dir  part1/main/tables_winsorise99 \
    --figures-dir part1/main/figures_winsorise99 \
    --cache-dir   part1/main/cache_winsorise99

echo "--- figures for winsorise99 ---"
conda run -n PhD python part1/main/manuscript/make_figures.py \
    --tables-dir part1/main/tables_winsorise99 \
    --out-dir    part1/main/manuscript/figures_winsorise99 \
    --cache-dir  part1/main/cache_winsorise99

# ------------------------------------------------------------------
# 5. Non-overlapping windows (stride = 3)
# ------------------------------------------------------------------
echo ""
echo "--- 5/5  --window-stride 3 ---"
conda run -n PhD python part1/main/main_analysis.py \
    --window-stride 3 \
    --tables-dir  part1/main/tables_stride3 \
    --figures-dir part1/main/figures_stride3 \
    --cache-dir   part1/main/cache_stride3

echo "--- figures for stride3 ---"
conda run -n PhD python part1/main/manuscript/make_figures.py \
    --tables-dir part1/main/tables_stride3 \
    --out-dir    part1/main/manuscript/figures_stride3 \
    --cache-dir  part1/main/cache_stride3

# ------------------------------------------------------------------
# Done
# ------------------------------------------------------------------
echo ""
echo "======================================================="
echo "  All sensitivities complete: $(date)"
echo "======================================================="
echo ""
echo "Results are in:"
echo "  part1/main/tables_health_board/"
echo "  part1/main/tables_size_offset/"
echo "  part1/main/tables_index_simd/"
echo "  part1/main/tables_winsorise99/"
echo "  part1/main/tables_stride3/"
echo ""
echo "Sensitivity figures are in:"
echo "  part1/main/manuscript/figures_health_board/"
echo "  part1/main/manuscript/figures_size_offset/"
echo "  part1/main/manuscript/figures_index_simd/"
echo "  part1/main/manuscript/figures_winsorise99/"
echo "  part1/main/manuscript/figures_stride3/"
