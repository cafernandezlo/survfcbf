#!/bin/bash
# Run one R script of the benchmark as a SLURM array in the CESGA apptainer container. The
# script path, relative to experiments/benchmark, is the first argument, and resources go
# on the sbatch command line. Each script's header gives the call used for the published
# run, e.g. from experiments/benchmark:
#   mkdir -p logs   # SLURM opens --output before the body runs
#   sbatch -p short --qos=short --time=03:00:00 --mem=32G --cpus-per-task=2 --array=1-11 \
#     cesga/submit_script.sbatch realdata/panel_redundancy.R
#SBATCH --job-name=survfcbf_script
#SBATCH --output=logs/script_%A_%a.out
#SBATCH --error=logs/script_%A_%a.out
set -euo pipefail
SCRIPT="${1:?give the R script, relative to experiments/benchmark}"

module load cesga/system
module load R/4.5.2-bioconductor-3.21-apptainer
export PATH="$PATH:/usr/bin:/bin"   # the apptainer wrapper needs coreutils on PATH
# the kit and its package library on the cluster: SURVFCBF_CLUSTER_ROOT, by default
# $STORE2/survfcbf on CESGA (STORE2 is the permanent per-user store there)
KIT="${SURVFCBF_CLUSTER_ROOT:-${STORE2:?set SURVFCBF_CLUSTER_ROOT to the directory holding experiments/benchmark and Rlib}/survfcbf}"
export R_LIBS_USER="${SURVFCBF_RLIB:-$KIT/Rlib}"
SIF=/opt/cesga/system/software/Core/R/4.5.2-bioconductor-3.21-apptainer/bioconductor_3.21-R-4.5.2.sif
APPT=/opt/cesga/system/software/Core/apptainer/1.2.3/bin/apptainer
BENCH="$KIT/experiments/benchmark"
export BENCH_DIR="$BENCH"           # the scripts setwd() here (apptainer cwd is unreliable)

cd "$BENCH"
export TMPDIR="${LUSTRE}/survfcbf_tmp/script_${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"
"$APPT" exec "$SIF" Rscript --no-init-file "$BENCH/$SCRIPT"
rc=$?
rm -rf "$TMPDIR"
exit $rc
