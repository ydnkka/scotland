#!/usr/bin/env bash
# run_chapter1.sh
# Run Chapter 1 primary analyses and core sensitivity analyses.
# Run from the repo (scotland/chapter1) with the PhD conda environment active:
#   conda activate PhD
#   bash run_chapter1.sh
# Or without activating first:
#   conda run -n PhD bash run_chapter1.sh

set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "======================================================="
echo "  Chapter 1 primary analyses"
echo "  $(date)"
echo "======================================================="

echo ""
echo "--- Main analyses ---"
conda run -n PhD python overall_analysis.py
echo "--- SIMD domain analyses ---"
conda run -n PhD python domain_analysis.py
echo "--- Wave analyses ---"
conda run -n PhD python wave_analysis.py


echo "======================================================="
echo "  Chapter 1 sensitivity analyses"
echo "======================================================="

# ------------------------------------------------------------------
# 1. Health-board spatial clustering
# ------------------------------------------------------------------
echo ""
echo "--- cluster-by health_board ---"
conda run -n PhD python overall_analysis.py \
    --cluster-by health_board \
    --tables-dir  sensitivity/tables_health_board \
    --figures-dir sensitivity/figures_health_board \
    --cache-dir   sensitivity/cache_health_board

# ------------------------------------------------------------------
# 2. Approximately non-overlapping windows
# ------------------------------------------------------------------
echo ""
echo "--- window-stride 3 ---"
conda run -n PhD python overall_analysis.py \
    --window-stride 3 \
    --tables-dir  sensitivity/tables_stride3 \
    --figures-dir sensitivity/figures_stride3 \
    --cache-dir   sensitivity/cache_stride3

# ------------------------------------------------------------------
# 3. Tail influence: winsorise top 1 percent
# ------------------------------------------------------------------
echo ""
echo "--- winsorise-quantile 0.99 ---"
conda run -n PhD python overall_analysis.py \
    --winsorise-quantile 0.99 \
    --tables-dir  sensitivity/tables_winsorise99 \
    --figures-dir sensitivity/figures_winsorise99 \
    --cache-dir   sensitivity/cache_winsorise99

# ------------------------------------------------------------------
# 4. Tail influence: exclude top 0.5 percent
# ------------------------------------------------------------------
echo ""
echo "--- exclude-tail-quantile 0.995 ---"
conda run -n PhD python overall_analysis.py \
    --exclude-tail-quantile 0.995 \
    --tables-dir  sensitivity/tables_exclude_tail995 \
    --figures-dir sensitivity/figures_exclude_tail995 \
    --cache-dir   sensitivity/cache_exclude_tail995

# ------------------------------------------------------------------
# Done
# ------------------------------------------------------------------
echo ""
echo "======================================================="
echo "  All complete: $(date)"
echo "======================================================="
echo ""
echo "Results are in:"
echo "  tables/"
echo "  sensitivity/tables_health_board/"
echo "  sensitivity/tables_stride3/"
echo "  sensitivity/tables_winsorise99/"
echo "  sensitivity/tables_exclude_tail995/"
