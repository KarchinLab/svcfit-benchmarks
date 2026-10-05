#!/usr/bin/env bash
#SBATCH --job-name=svcfit_3arm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-2249%25
set -euo pipefail
: "${THREE_ARM_MANIFEST:?}" "${THREE_ARM_OUTPUT_ROOT:?}" "${VISOR_PIPELINE_DIR:?}"
: "${WORKFLOW_REPO:?}" "${EXPECTED_WORKFLOW_COMMIT:?}" "${EXPECTED_SVCFIT_COMMIT:?}"
: "${SVCFIT_RSCRIPT:?}" "${SVCFIT_PKG_DIR:?}" "${SVCFIT_TRUTH_DIR:?}" "${SVCFIT_RLIB:?}"
: "${SLURM_ARRAY_TASK_ID:?}"
[[ "$(git -C "$WORKFLOW_REPO" rev-parse HEAD)" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
  echo "ERROR: workflow checkout moved after submission" >&2; exit 1;
}
[[ "$(git -C "$SVCFIT_PKG_DIR" rev-parse HEAD)" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
  echo "ERROR: SVCFit checkout moved after submission" >&2; exit 1;
}
"$SVCFIT_RSCRIPT" "$VISOR_PIPELINE_DIR/08_extract_svcfit_condition.R" \
  --manifest "$THREE_ARM_MANIFEST" --task-id "$SLURM_ARRAY_TASK_ID" \
  --output-root "$THREE_ARM_OUTPUT_ROOT" --svcfit-package "$SVCFIT_PKG_DIR" \
  --truth-dir "$SVCFIT_TRUTH_DIR" --rlib "$SVCFIT_RLIB"
