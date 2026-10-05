#!/bin/bash
#SBATCH --job-name=a2_extract
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/a2_extract_%A_%a.out
#SBATCH --error=logs/a2_extract_%A_%a.err
#SBATCH --array=0-79
#
# 4 purities x 20 reps = 80 array tasks.
# Each task processes both pre and post timepoints sequentially in one R startup.
# Edit --account, --partition, --time, --mem, and the `module load` line below.

set -euo pipefail

module load r/4.3.0   # FILL IN -- your cluster's R module

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export A2_SCRIPT_DIR="$SCRIPT_DIR"

mkdir -p logs

PURITIES=(20 40 60 80)
REPS=($(seq 1 20))

idx=$SLURM_ARRAY_TASK_ID
n_rep=${#REPS[@]}

rep_idx=$(( idx % n_rep ))
pur_idx=$(( idx / n_rep ))

PURITY=${PURITIES[$pur_idx]}
REP=${REPS[$rep_idx]}

echo "[$(date)] Task $idx: purity=$PURITY rep=$REP (both timepoints)"

for TP in pre post; do
  echo "  -> $TP"
  Rscript "$SCRIPT_DIR/extract_per_run.R" "$PURITY" "$REP" "$TP"
done

echo "[$(date)] Task $idx complete."
