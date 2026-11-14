#!/bin/bash
# One grid cell per array task: the array INDEX is the cell index of grid_cells(), because
# N_TASKS is left unset and run_cell.R then maps task i to cell i. This lets cells be
# grouped by expected runtime and sent to different QoS (cesga_core cells range from
# ~6 min at p=200 to ~5 h at n=1200 on epyc nodes, and the *_short QoS cap walltime at
# 6 h). Cells are single-threaded (TotalCPU = Elapsed in the August 2026 run) and peak
# below 3 GB, hence 1 CPU / 4 GB. In cesga_core the cells of scenario s (1-11, in
# grid_cells() order) are s, s+11, ..., s+11*29. The November 2026 rerun used three
# groups, from experiments/benchmark:
#   mkdir -p logs results/cesga_core/cells   # SLURM opens --output before the body runs
#   sbatch -p short --qos=short              --time=02:00:00 --array=<n150,n300,p200 cells> cesga/submit_grid_cell.sbatch cesga_core
#   sbatch -p short --qos=clk_short  -C clk  --time=06:00:00 --array=<other cells>         cesga/submit_grid_cell.sbatch cesga_core
#   sbatch -p short --qos=epyc_short -C epyc --time=06:00:00 --array=<n1200 cells>         cesga/submit_grid_cell.sbatch cesga_core
# Read the predicted start with `sbatch --test-only` before choosing: on 2026-11-04 any job
# over 6 h in medium was predicted six days out, while the 6 h *_short groups started the
# same day. Use submit_grid.sbatch instead for the strided multi-cell layout.
#SBATCH --job-name=survfcbf_cell
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
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
unset N_TASKS               # one cell per task: run_cell.R maps array index i to cell i

cd "$BENCH"
mkdir -p logs "results/${TIER}/cells"

# scratch workdir per task in LUSTRE (regenerated each run, never the permanent store)
export TMPDIR="${LUSTRE}/survfcbf_tmp/${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"

rc=0   # under set -e a bare failure would exit here and leave the scratch dir behind
"$APPT" exec "$SIF" Rscript --no-init-file "$BENCH/cesga/run_cell.R" "$TIER" || rc=$?
rm -rf "$TMPDIR"
exit $rc
