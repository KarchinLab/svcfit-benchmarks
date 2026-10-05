#!/bin/bash
#SBATCH --job-name=a2_agg
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/a2_agg_%j.out
#SBATCH --error=logs/a2_agg_%j.err
#
# Single-job aggregation + statistics + Fig A2-3.
# Submit with afterok dependency on the extract array (see run_all.sh).

set -euo pipefail

module load r/4.3.0   # FILL IN -- your cluster's R module

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export A2_SCRIPT_DIR="$SCRIPT_DIR"

mkdir -p logs

echo "[$(date)] Aggregating per-run outputs and running statistics."
Rscript "$SCRIPT_DIR/aggregate_and_analyze.R"
echo "[$(date)] Done."
