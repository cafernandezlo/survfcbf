#!/bin/bash
# CESGA SLURM array for the simulation benchmark grid. Each array
# task runs a STRIDED subset of grid cells via cesga/run_cell.R inside the
# bioconductor R 4.5.2 apptainer container; run_cell.R reads N_TASKS (exported
# here from SLURM_ARRAY_TASK_COUNT) to stride. Keep the array <= 100 tasks: the
# short QoS caps queued jobs per user at MaxSubmitPU=100, so one task per cell
# (330 cells) is rejected. ~80 tasks x ~4 cells is a good split. Submit, e.g.:
#   mkdir -p logs results/cesga_core/cells   # SLURM opens --output before the body runs
#   sbatch --array=1-80%40 cesga/submit_grid.sbatch cesga_core
# Run from experiments/benchmark. Aggregate afterwards with cesga/aggregate.R.
#SBATCH --job-name=survfcbf_grid
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/cell_%A_%a.out
#SBATCH --error=logs/cell_%A_%a.out
set -euo pipefail

TIER="${1:-cesga_core}"

module load cesga/system
module load R/4.5.2-bioconductor-3.21-apptainer
export PATH="$PATH:/usr/bin:/bin"   # the apptainer wrapper needs coreutils on PATH (module loads can drop them)
# the kit and its package library on the cluster: SURVFCBF_CLUSTER_ROOT, by default
# $STORE2/survfcbf on CESGA (STORE2 is the permanent per-user store there)
KIT="${SURVFCBF_CLUSTER_ROOT:-${STORE2:?set SURVFCBF_CLUSTER_ROOT to the directory holding experiments/benchmark and Rlib}/survfcbf}"
export R_LIBS_USER="${SURVFCBF_RLIB:-$KIT/Rlib}"
SIF=/opt/cesga/system/software/Core/R/4.5.2-bioconductor-3.21-apptainer/bioconductor_3.21-R-4.5.2.sif
APPT=/opt/cesga/system/software/Core/apptainer/1.2.3/bin/apptainer
BENCH="$KIT/experiments/benchmark"
export BENCH_DIR="$BENCH"   # run_cell.R setwd()s here (apptainer container cwd is unreliable)
# total array size -> run_cell.R uses it to stride cells across tasks (unset for non-array runs)
if [ -n "${SLURM_ARRAY_TASK_COUNT:-}" ]; then export N_TASKS="$SLURM_ARRAY_TASK_COUNT"; fi

cd "$BENCH"
mkdir -p logs "results/${TIER}/cells"

# scratch workdir per task in LUSTRE (regenerated each run; never the permanent store)
export TMPDIR="${LUSTRE}/survfcbf_tmp/${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"

# absolute script path (container cwd is unreliable); run_cell.R setwd()s to BENCH_DIR
"$APPT" exec "$SIF" Rscript --no-init-file "$BENCH/cesga/run_cell.R" "$TIER"
rc=$?
rm -rf "$TMPDIR"
exit $rc
