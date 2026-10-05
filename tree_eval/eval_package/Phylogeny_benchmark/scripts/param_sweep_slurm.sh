#!/bin/bash
#SBATCH --job-name=param_sweep
#SBATCH --output=log/param_sweep_%A_%a.out
#SBATCH --error=log/param_sweep_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --time=1:00:00
#SBATCH --array=0-179

# 10 test cases × 6 concentration values × 3 min_dist = 180 tasks
# Pipeline: svcfit → clustering → tree (all output inside param_sweep/)
# Submit: sbatch param_sweep_slurm.sh
# Merge:  Rscript param_sweep.R --merge

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
# shellcheck disable=SC1091
source "${VISOR_CONFIG:-$REPO_ROOT/config.local.sh}"
PUB_DIR="$(dirname "$SCRIPT_DIR")"
SCRP_DIR=$SCRIPT_DIR
SWEEP_DIR=$PUB_DIR/output/param_sweep

# Task 0 removes previous results; all others wait for it to finish
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    rm -rf "$SWEEP_DIR"
    mkdir -p "$SWEEP_DIR"
else
    sleep 10
fi

"$RSCRIPT" "$SCRP_DIR/param_sweep.R" --task "$SLURM_ARRAY_TASK_ID"
