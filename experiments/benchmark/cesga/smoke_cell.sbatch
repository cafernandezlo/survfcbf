#!/bin/bash
# One-cell smoke for the corrected-recovery re-run. Forces N_TASKS to the full cell
# count so the strided driver picks exactly one cell (SMOKE_CELL), instead of the
# whole tail that --array=k-k would give (SLURM_ARRAY_TASK_COUNT would be 1).
#SBATCH --job-name=survfcbf_smoke
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/smoke_%j.out
#SBATCH --error=logs/smoke_%j.out
set -euo pipefail

TIER="${1:-cesga_core}"
CELL="${2:-9}"

module load cesga/system
module load R/4.5.2-bioconductor-3.21-apptainer
export PATH="$PATH:/usr/bin:/bin"
# the kit and its package library on the cluster: SURVFCBF_CLUSTER_ROOT, by default
# $STORE2/survfcbf on CESGA (STORE2 is the permanent per-user store there)
KIT="${SURVFCBF_CLUSTER_ROOT:-${STORE2:?set SURVFCBF_CLUSTER_ROOT to the directory holding experiments/benchmark and Rlib}/survfcbf}"
export R_LIBS_USER="${SURVFCBF_RLIB:-$KIT/Rlib}"
SIF=/opt/cesga/system/software/Core/R/4.5.2-bioconductor-3.21-apptainer/bioconductor_3.21-R-4.5.2.sif
APPT=/opt/cesga/system/software/Core/apptainer/1.2.3/bin/apptainer
BENCH="$KIT/experiments/benchmark"
export BENCH_DIR="$BENCH"

# 330 cells total; id=CELL with stride 330 selects exactly cell CELL
export N_TASKS=330
export SLURM_ARRAY_TASK_ID="$CELL"

cd "$BENCH"
mkdir -p logs "results/${TIER}/cells"

export TMPDIR="${LUSTRE}/survfcbf_tmp/smoke_${SLURM_JOB_ID:-0}"
mkdir -p "$TMPDIR"

"$APPT" exec "$SIF" Rscript --no-init-file "$BENCH/cesga/run_cell.R" "$TIER"
rc=$?
rm -rf "$TMPDIR"
exit $rc
