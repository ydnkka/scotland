#!/usr/bin/env bash
# run_chapter1.sh
# Run all five Part 1 sensitivity analyses and regenerate figures.
# Run from the repo (scotland/part1) with the PhD conda environment active:
#   conda activate PhD
#   bash run_chapter1.sh
# Or without activating first:
#   conda run -n PhD bash run_chapter1.sh

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


echo "======================================================="
echo "  Part 1 sensitivity analyses"
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
# Done
# ------------------------------------------------------------------
echo ""
echo "======================================================="
echo "  All Complete: $(date)"
echo "======================================================="
echo ""
echo "Results are in:"
echo "  tables/"
echo "  sensitivity/tables_health_board/"


#echo ""
#echo "======================================================="
#echo "  Make figures for all analyses (primary + sensitivity)"
#echo "======================================================="

#echo "--- Primary analysis figures ---"
#conda run -n PhD python manuscript/make_figures.py

#echo "--- figures for health_board ---"
#conda run -n PhD python manuscript/make_figures.py \
#    --tables-dir sensitivity/tables_health_board \
#    --out-dir    manuscript/sensitivity/figures_health_board \
#    --cache-dir  sensitivity/cache_health_board

#echo "--- figures for winsorise99 ---"
#conda run -n PhD python manuscript/make_figures.py \
#    --tables-dir sensitivity/tables_winsorise99 \
#    --out-dir    manuscript/sensitivity/figures_winsorise99 \
#    --cache-dir  sensitivity/cache_winsorise99


#echo ""
#echo "Figures are in:"
#echo "  manuscript/figures/"
#echo "  manuscript/sensitivity/figures_health_board/"
#echo "  manuscript/sensitivity/figures_winsorise99/"
