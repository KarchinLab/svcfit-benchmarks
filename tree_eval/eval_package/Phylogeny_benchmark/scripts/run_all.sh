#!/bin/bash
###############################################################################
# run_all.sh
#
# Master script to run the full longitudinal simulation benchmark.
# Submits SLURM jobs in dependency order for all scenario/purity combinations
# across five simulation replicates (BOOT=0..4), each with its own seed range
# (see longi_short.sh): VISOR SHORtS -> Manta/SVtyper/FACETS -> SVCFit.
#
# Replicates per scenario/purity: 5 replicates x 5 SV-CNV configurations = 25.
#
# The reported S1 results re-run SVCFit, clustering and tree reconstruction at
# SVCFit 7f32d81 on these simulations with run_svcfit_and_evaluate.sh --scope full.
#
# Usage:
#   bash run_all.sh
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
# shellcheck disable=SC1091
source "${VISOR_CONFIG:-$REPO_ROOT/config.local.sh}"
PUB_DIR="$(dirname "$SCRIPT_DIR")"
SCENARIOS=(S1 S2 S3 S4)
PURITIES=(10 20 40 60 80)
BOOTS=(0 1 2 3 4)

# --- PATHS ---
WORK_DIR=$PUB_DIR/outputs
TRUTH_DIR=$PUB_DIR/data/hack
EVAL_DIR=$PUB_DIR/outputs/evaluation

mkdir -p log

echo "=== Longitudinal Simulation Benchmark (5 simulation replicates) ==="
echo "Scenarios: ${SCENARIOS[@]}"
echo "Purities:  ${PURITIES[@]}"
echo "Boots:     ${BOOTS[@]}"
echo ""

LAST_JOBS=""

for BOOT in "${BOOTS[@]}"; do
  echo "=== Simulation replicate BOOT=$BOOT ==="

  for SCENARIO in "${SCENARIOS[@]}"; do
    for PPUR in "${PURITIES[@]}"; do
      echo "  --- Submitting: $SCENARIO / purity=$PPUR / boot=$BOOT ---"

        JOB1=$(sbatch --parsable \
          --export=ALL,SCENARIO=$SCENARIO,ppur=$PPUR,BOOT=$BOOT \
          $SCRIPT_DIR/longi_short.sh)
        echo "  Step 1 (VISOR SHORtS): job $JOB1"

        JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 \
          --export=ALL,SCENARIO=$SCENARIO,ppur=$PPUR,BOOT=$BOOT \
          $SCRIPT_DIR/longi_calling.sh)
        echo "  Step 2 (pipeline):     job $JOB2"

        JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 \
          --export=ALL,SCENARIO=$SCENARIO,ppur=$PPUR,BOOT=$BOOT \
          $SCRIPT_DIR/longi_svcfit.sh)
        echo "  Step 3 (SVCFit):       job $JOB3"
        LAST_JOBS="$LAST_JOBS:$JOB3"

    done
  done
done

# Evaluation runs after ALL scenario/purity/boot jobs complete.
LAST_JOBS=${LAST_JOBS#:}  # strip leading colon
echo ""
echo "--- Submitting evaluation (depends on all above) ---"
JOB_EVAL=$(sbatch --parsable --dependency=afterok:$LAST_JOBS \
  --job-name=longi_eval \
  --output=log/longi_eval_%j.out \
  --error=log/longi_eval_%j.err \
  --nodes=1 --mem=4G --time=15:00 \
  --wrap="'$RSCRIPT' '$SCRIPT_DIR/evaluate_downstream.R' \
      --work_dir $WORK_DIR \
      --truth_dir $TRUTH_DIR \
      --out_dir $EVAL_DIR")
echo "  Evaluation: job $JOB_EVAL"

echo ""
echo "=== All jobs submitted ==="
echo "Monitor with: squeue -u \$USER"
echo "Evaluation output will be in: $EVAL_DIR"
