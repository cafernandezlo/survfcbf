#!/bin/bash
# Generic one-cohort-per-task SLURM array for the real-data stage. Every
# array script here (run_dataset, extract_signatures, extract_scores, functional_stability,
# cardinality_curve, car_ablation, module_*) reads SLURM_ARRAY_TASK_ID, indexes the cohort
# list and writes one output per cohort with skip-if-exists, so they all share this
# launcher instead of a per-script copy. Time and memory come from the environment, so a
# heavy stage does not force everyone else onto a big allocation.
#
# Usage (from experiments/benchmark, after realdata/Data is populated):
#   mkdir -p realdata/logs
#   RD_SCRIPT=realdata/functional_stability.R sbatch --array=1-11 realdata/submit_array.sbatch
#   RD_SCRIPT=realdata/cardinality_curve.R RD_TIME=04:00:00 RD_MEM=32G \
#     sbatch --array=1-11 realdata/submit_array.sbatch
#
# Defaults suit the light probes (1 h, 32 G). The full benchmark (run_dataset.R) needs
# RD_TIME=20:00:00 RD_MEM=64G -- see submit_realdata.sbatch, kept separate because its
# ~12-14 h per cohort is an overnight run, not a probe.
#SBATCH --job-name=survfcbf_array
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=realdata/logs/%x_%A_%a.out
#SBATCH --error=realdata/logs/%x_%A_%a.out
set -euo pipefail

: "${RD_SCRIPT:?set RD_SCRIPT to the array script, e.g. realdata/functional_stability.R}"

module load cesga/system
module load R/4.5.2-bioconductor-3.21-apptainer
export PATH="$PATH:/usr/bin:/bin"          # the apptainer wrapper needs coreutils on PATH
# the kit and its package library on the cluster: SURVFCBF_CLUSTER_ROOT, by default
# $STORE2/survfcbf on CESGA (STORE2 is the permanent per-user store there)
KIT="${SURVFCBF_CLUSTER_ROOT:-${STORE2:?set SURVFCBF_CLUSTER_ROOT to the directory holding experiments/benchmark and Rlib}/survfcbf}"
export R_LIBS_USER="${SURVFCBF_RLIB:-$KIT/Rlib}"
SIF=/opt/cesga/system/software/Core/R/4.5.2-bioconductor-3.21-apptainer/bioconductor_3.21-R-4.5.2.sif
APPT=/opt/cesga/system/software/Core/apptainer/1.2.3/bin/apptainer
BENCH="$KIT/experiments/benchmark"
export BENCH_DIR="$BENCH"                  # scripts setwd() here (apptainer cwd is unreliable)
export BOMMERT_DATA="$BENCH/realdata/Data" # loaders read the cached ARFFs offline
export METABRIC_DIR="$BENCH/realdata/METABRIC"

cd "$BENCH"
mkdir -p realdata/logs
export TMPDIR="${LUSTRE}/survfcbf_tmp/${SLURM_JOB_NAME:-arr}_${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"

"$APPT" exec "$SIF" Rscript --no-init-file "$BENCH/$RD_SCRIPT"
rc=$?
rm -rf "$TMPDIR"
exit $rc
