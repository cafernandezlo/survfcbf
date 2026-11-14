#!/bin/bash
# Single entry point for the CESGA half of the real-data stage. Submits
# every one-cohort-per-task array through the generic launcher, with the time/memory each
# stage actually needs. All the array scripts write one output per cohort with
# skip-if-exists, so re-running this only fills the gaps.
#
# Usage (from experiments/benchmark on a CESGA login node, after realdata/Data is cached):
#   bash realdata/submit_all_cesga.sh            # everything except the overnight benchmark
#   bash realdata/submit_all_cesga.sh --with-benchmark
#   bash realdata/submit_all_cesga.sh --dry-run
#
# The main benchmark (run_dataset.R) is opt-in: ~12-14 h per cohort, an overnight run, and
# it uses its own submit_realdata.sbatch (medium QoS, 20 h, 64 G).
#
# Scheduling note: these arrays compete for nodes via FairShare. If another project of the
# same user has many jobs RUNNING, these sit in PENDING(Priority). Freeze that project's
# PENDING array with `scontrol hold <jobid>` (reversible with `release`) and release it as
# soon as these flip to RUNNING -- once running they hold their nodes, so keeping the other
# project held any longer only starves it.
set -euo pipefail
cd "$(dirname "$0")/.."                     # experiments/benchmark

DRY=0; WITH_BENCH=0
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    --with-benchmark) WITH_BENCH=1 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done
mkdir -p realdata/logs

submit () {                                  # submit <script> <jobname> <time> <mem>
  local script="$1" name="$2" time="$3" mem="$4"
  if [ ! -f "$script" ]; then echo "MISSING $script" >&2; return 1; fi
  local cmd=(sbatch --parsable --array=1-11 --job-name="$name" --time="$time" --mem="$mem"
             --export="ALL,RD_SCRIPT=$script" realdata/submit_array.sbatch)
  if [ "$DRY" = "1" ]; then printf 'DRY  %-28s %s %s\n' "$name" "$time" "$mem"; return 0; fi
  local jid; jid=$("${cmd[@]}")
  printf 'submitted %-28s job %s (%s, %s)\n' "$name" "$jid" "$time" "$mem"
}

# Stage 1: signatures and per-gene scores (feed the enrichment/GSEA analyses)
submit realdata/extract_signatures.R  sf_signatures 01:00:00 32G
CODING_ONLY=1 submit realdata/extract_signatures.R sf_signatures_coding 01:00:00 32G
submit realdata/extract_scores.R      sf_scores     01:00:00 32G

# Stage 2: the evidence-gap tests
submit realdata/functional_stability.R sf_functional 01:00:00 32G
submit realdata/cardinality_curve.R    sf_cardinality 04:00:00 32G

# Stage 3: mechanism and the module extension
submit realdata/car_ablation.R        sf_car        01:00:00 32G
submit realdata/module_probe_v3.R     sf_modules    01:00:00 32G
submit realdata/module_predict.R      sf_mod_pred   01:00:00 32G
submit realdata/module_leaksafe.R     sf_mod_leak   02:00:00 32G

# Stage 4 (opt-in): the full nested-CV benchmark, ~12-14 h per cohort
if [ "$WITH_BENCH" = "1" ]; then
  if [ "$DRY" = "1" ]; then echo "DRY  benchmark (submit_realdata.sbatch, 20 h, 64 G)"
  else
    jid=$(sbatch --parsable --array=1-11 realdata/submit_realdata.sbatch)
    echo "submitted benchmark                  job $jid (20:00:00, 64G)"
  fi
else
  echo "skipped   benchmark (run_dataset.R)   -- pass --with-benchmark for the overnight run"
fi

cat <<'EOF'

Next:
  squeue -u "$USER"                       # watch
  # when the arrays finish, copy results/ back and run the local half:
  #   Rscript --no-init-file realdata/run_all_local.R
EOF
