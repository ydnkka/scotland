#!/usr/bin/env bash
# part1_sensitivities.sh
# Run all five Part 1 sensitivity analyses and regenerate figures.
# Run from the repo (scotland/part1) with the PhD conda environment active:
#   conda activate PhD
#   bash part1_sensitivities.sh
# Or without activating first:
#   conda run -n PhD bash part1_sensitivities.sh

set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "======================================================="
echo "  Part 1 primary analyses"
echo "  $(date)"
echo "======================================================="

echo ""
echo "--- Main analyses ---"
conda run -n PhD python overall_analysis.py
echo "--- SIMD domain analyses ---"
conda run -n PhD python domain_analysis.py
echo "--- Wave analyses ---"
conda run -n PhD python wave_analysis.py
echo "--- Log-linear sensitivity analyses ---"
conda run -n PhD python loglinear_sensitivity.py
echo "--- Primary analysis figures ---"
conda run -n PhD python manuscript/make_figures.py


echo "======================================================="
echo "  Part 1 sensitivity analyses"
echo "======================================================="

# ------------------------------------------------------------------
# 1. Health-board spatial clustering
# ------------------------------------------------------------------
echo ""
echo "--- 1/5  --cluster-by health_board ---"
conda run -n PhD python overall_analysis.py \
    --cluster-by health_board \
    --tables-dir  sensitivity/tables_health_board \
    --figures-dir sensitivity/figures_health_board \
    --cache-dir   sensitivity/cache_health_board

echo "--- figures for health_board ---"
conda run -n PhD python manuscript/make_figures.py \
    --tables-dir sensitivity/tables_health_board \
    --out-dir    manuscript/sensitivity/figures_health_board \
    --cache-dir  sensitivity/cache_health_board

# ------------------------------------------------------------------
# 2. Sampling-pool size offset (wn_no_sequences)
# ------------------------------------------------------------------
echo ""
echo "--- 2/5  --use-size-offset ---"
conda run -n PhD python overall_analysis.py \
    --use-size-offset \
    --tables-dir  sensitivity/tables_size_offset \
    --figures-dir sensitivity/figures_size_offset \
    --cache-dir   sensitivity/cache_size_offset

echo "--- figures for size_offset ---"
conda run -n PhD python manuscript/make_figures.py \
    --tables-dir sensitivity/tables_size_offset \
    --out-dir    manuscript/sensitivity/figures_size_offset \
    --cache-dir  sensitivity/cache_size_offset

# ------------------------------------------------------------------
# 3. Index-case SIMD instead of mean-cluster SIMD
# ------------------------------------------------------------------
echo ""
echo "--- 3/5  --use-index-simd ---"
conda run -n PhD python overall_analysis.py \
    --use-index-simd \
    --tables-dir  sensitivity/tables_index_simd \
    --figures-dir sensitivity/figures_index_simd \
    --cache-dir   sensitivity/cache_index_simd

echo "--- figures for index_simd ---"
conda run -n PhD python manuscript/make_figures.py \
    --tables-dir sensitivity/tables_index_simd \
    --out-dir    manuscript/sensitivity/figures_index_simd \
    --cache-dir  sensitivity/cache_index_simd

# ------------------------------------------------------------------
# 4. Tail winsorisation at 99th percentile
# ------------------------------------------------------------------
echo ""
echo "--- 4/5  --winsorise-quantile 0.99 ---"
conda run -n PhD python overall_analysis.py \
    --winsorise-quantile 0.99 \
    --tables-dir  sensitivity/tables_winsorise99 \
    --figures-dir sensitivity/figures_winsorise99 \
    --cache-dir   sensitivity/cache_winsorise99

echo "--- figures for winsorise99 ---"
conda run -n PhD python manuscript/make_figures.py \
    --tables-dir sensitivity/tables_winsorise99 \
    --out-dir    manuscript/sensitivity/figures_winsorise99 \
    --cache-dir  sensitivity/cache_winsorise99

# ------------------------------------------------------------------
# 5. Non-overlapping windows (stride = 3)
# ------------------------------------------------------------------
echo ""
echo "--- 5/5  --window-stride 3 ---"
conda run -n PhD python overall_analysis.py \
    --window-stride 3 \
    --tables-dir  sensitivity/tables_stride3 \
    --figures-dir sensitivity/figures_stride3 \
    --cache-dir   sensitivity/cache_stride3

echo "--- figures for stride3 ---"
conda run -n PhD python manuscript/make_figures.py \
    --tables-dir sensitivity/tables_stride3 \
    --out-dir    manuscript/sensitivity/figures_stride3 \
    --cache-dir  sensitivity/cache_stride3

# ------------------------------------------------------------------
# Done
# ------------------------------------------------------------------
echo ""
echo "======================================================="
echo "  All sensitivities complete: $(date)"
echo "======================================================="
echo ""
echo "Results are in:"
echo "  sensitivity/tables_health_board/"
echo "  sensitivity/tables_size_offset/"
echo "  sensitivity/tables_index_simd/"
echo "  sensitivity/tables_winsorise99/"
echo "  sensitivity/tables_stride3/"
echo ""
echo "Sensitivity figures are in:"
echo "  manuscript/sensitivity/figures_health_board/"
echo "  manuscript/sensitivity/figures_size_offset/"
echo "  manuscript/sensitivity/figures_index_simd/"
echo "  manuscript/sensitivity/figures_winsorise99/"
echo "  manuscript/sensitivity/figures_stride3/"
