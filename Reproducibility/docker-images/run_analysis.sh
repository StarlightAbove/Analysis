#!/usr/bin/env bash
#
# run_analysis.sh — ENTRYPOINT for cnv-analysis (Stage 2).
#
# Runs the three R scripts headless and in order. `set -e` halts the chain the
# moment any stage fails, so run_all.R never runs on incomplete tool output.
set -euo pipefail

cd /Analysis

echo "[1/3] Conumee arm ..."
Rscript /Analysis/Scripts/conumee.R

echo "[2/3] Sesame arm ..."
Rscript /Analysis/Scripts/sesame.R

echo "[3/3] Analysis (Analytic Code/run_all.R) ..."
Rscript "/Analysis/Analytic Code/run_all.R"

echo "Done. Results in /Analysis/Outputs and /Analysis/results."
