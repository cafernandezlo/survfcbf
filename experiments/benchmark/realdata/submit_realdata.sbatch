#!/bin/bash
# SLURM array for the Bommert real-data benchmark. One task per
# cohort (11), each loads its cached ARFFs and runs the harness. p~23k makes each
# cell heavy, so this uses the medium QoS (the 20h --time routes there; 11 tasks is
# well under medium's MaxSubmitPU=50) instead of the 6h short cap. Submit from
# experiments/benchmark, AFTER the data is cached (realdata/download_all.R on login):
#   mkdir -p realdata/logs results/realdata/cells
#   sbatch --array=1-11 realdata/submit_realdata.sbatch
# Aggregate afterwards with realdata/aggregate_realdata.R.
#SBATCH --job-name=survfcbf_realdata
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=realdata/logs/ds_%A_%a.out
#SBATCH --error=realdata/logs/ds_%A_%a.out
set -euo pipefail

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
export BENCH_DIR="$BENCH"                 # run_dataset.R setwd()s here (apptainer cwd is unreliable)
export BOMMERT_DATA="$BENCH/realdata/Data" # loader reads the pre-downloaded cache offline

cd "$BENCH"
mkdir -p realdata/logs results/realdata/cells

export TMPDIR="${LUSTRE}/survfcbf_tmp/rd_${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"

"$APPT" exec "$SIF" Rscript --no-init-file "$BENCH/realdata/run_dataset.R"
rc=$?
rm -rf "$TMPDIR"
exit $rc
