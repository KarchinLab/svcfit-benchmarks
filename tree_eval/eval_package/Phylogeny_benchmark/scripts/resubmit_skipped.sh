#!/bin/bash
###############################################################################
# resubmit_skipped.sh
#
# Scans all boot/scenario/purity/exp combinations and resubmits only the
# cases that are missing from the evaluation (i.e., skipped by
# evaluate_downstream.R due to absent output files).
#
# Assigns the minimum required stages based on what already exists on disk:
#   tree               — clustering done, only tree/tree_result.rds missing
#   cluster+tree       — svcfit_output exists, clustering/tree both missing
#   svcfit+cluster+tree — only calling output (t1/t2) present
#
# Groups missing experiments by (boot, scenario, purity, stages) and submits
# one sbatch array job per group so that overhead is minimised.
#
# Usage:
#   bash resubmit_skipped.sh [--dry-run]
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
# shellcheck disable=SC1091
source "${VISOR_CONFIG:-$REPO_ROOT/config.local.sh}"
PUB_DIR="$(dirname "$SCRIPT_DIR")"
WORK_DIR=$PUB_DIR/outputs
TRUTH_DIR=$PUB_DIR/data/hack
EVAL_DIR=$WORK_DIR/evaluation/boots0_2
COV=50

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

mkdir -p "$SCRIPT_DIR/log"

scenarios=(S1 S2 S3 S4)
purities=(10 20 40 60 80)
boots=(0 1 2)
exp_lst=(exp1 exp2 exp3 exp4 exp5)

echo "=== Scanning for skipped cases (COV=$COV) ==="
[[ $DRY_RUN -eq 1 ]] && echo "    (dry-run mode: no jobs submitted)"
echo ""

LAST_JOBS=""
total_submitted=0

for boot in "${boots[@]}"; do
  for scen in "${scenarios[@]}"; do
    for pur in "${purities[@]}"; do

      # Collect missing array indices per stages tier
      indices_tree=""
      indices_cluster=""
      indices_svcfit=""

      for i in 0 1 2 3 4; do
        exp=${exp_lst[$i]}

        if [[ $boot -eq 0 ]]; then
          dir="$WORK_DIR/$scen/c${COV}p${pur}/$exp"
        else
          dir="$WORK_DIR/$scen/c${COV}p${pur}/b${boot}/$exp"
        fi

        cluster_ok=0; tree_ok=0; svcfit_ok=0
        [[ -f "$dir/clustering/sv2cluster.csv" && \
           -f "$dir/clustering/cluster_centroids.csv" ]] && cluster_ok=1
        [[ -f "$dir/tree/tree_result.rds" ]] && tree_ok=1
        [[ -d "$dir/svcfit_output" ]] && svcfit_ok=1

        [[ $tree_ok -eq 1 ]] && continue   # nothing to do

        if [[ $cluster_ok -eq 1 ]]; then
          indices_tree="${indices_tree:+$indices_tree,}$i"
        elif [[ $svcfit_ok -eq 1 ]]; then
          indices_cluster="${indices_cluster:+$indices_cluster,}$i"
        else
          indices_svcfit="${indices_svcfit:+$indices_svcfit,}$i"
        fi
      done

      # Submit one job per non-empty tier
      for tier in "tree:$indices_tree" "cluster+tree:$indices_cluster" "svcfit+cluster+tree:$indices_svcfit"; do
        stages="${tier%%:*}"
        indices="${tier#*:}"
        [[ -z "$indices" ]] && continue

        echo "  boot=$boot $scen/p${pur} STAGES=$stages array=$indices"
        total_submitted=$((total_submitted + 1))

        if [[ $DRY_RUN -eq 0 ]]; then
          JOB=$(sbatch --parsable \
            --array="$indices" \
            --export=ALL,SCENARIO=$scen,ppur=$pur,BOOT=$boot,STAGES=$stages \
            "$SCRIPT_DIR/longi_svcfit.sh")
          echo "    -> job $JOB"
          LAST_JOBS="${LAST_JOBS:+$LAST_JOBS:}$JOB"
        fi
      done

    done
  done
done

echo ""
echo "=== ${total_submitted} job(s) $([ $DRY_RUN -eq 1 ] && echo 'would be submitted' || echo 'submitted') ==="

if [[ $DRY_RUN -eq 0 && -n "$LAST_JOBS" ]]; then
  echo ""
  echo "--- Submitting evaluation (depends on all above) ---"
  JOB_EVAL=$(sbatch --parsable \
    --dependency=afterok:$LAST_JOBS \
    --job-name=longi_eval_rerun \
    --output=log/longi_eval_rerun_%j.out \
    --error=log/longi_eval_rerun_%j.err \
    --nodes=1 --mem=4G --time=15:00 \
    --wrap="'$RSCRIPT' '$SCRIPT_DIR/evaluate_downstream.R' \
        --work_dir $WORK_DIR \
        --truth_dir $TRUTH_DIR \
        --out_dir $EVAL_DIR")
  echo "  Evaluation: job $JOB_EVAL"
fi

echo ""
echo "Monitor with: squeue -u \$USER"
