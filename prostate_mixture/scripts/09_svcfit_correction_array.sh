#!/usr/bin/env bash
#SBATCH --job-name=pm_svcf_shadow
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-29%10
set -euo pipefail

: "${PROSTATE_OUTPUT_ROOT:?}" "${PROSTATE_DATA_ROOT:?}" "${PROSTATE_SCRIPT_DIR:?}"
: "${WORKFLOW_REPO:?}" "${EXPECTED_WORKFLOW_COMMIT:?}"
: "${SVCFIT_PKG_DIR:?}" "${EXPECTED_SVCFIT_COMMIT:?}"
: "${SVCFIT_RSCRIPT:?}" "${SVCFIT_RLIB:?}" "${SVCFIT_RLIB_COMMIT_FILE:?}"
: "${SLURM_ARRAY_TASK_ID:?}"

[[ "$(git -C "$WORKFLOW_REPO" rev-parse HEAD)" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
  echo "ERROR: workflow checkout moved after submission" >&2; exit 1;
}
[[ "$(git -C "$SVCFIT_PKG_DIR" rev-parse HEAD)" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
  echo "ERROR: SVCFit checkout moved after submission" >&2; exit 1;
}
[[ "$(tr -d '[:space:]' < "$SVCFIT_RLIB_COMMIT_FILE")" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
  echo "ERROR: installed SVCFit library provenance does not match the pinned commit" >&2; exit 1;
}

rep=$((SLURM_ARRAY_TASK_ID + 1))
parts="$PROSTATE_OUTPUT_ROOT/parts"
mkdir -p "$parts" "$PROSTATE_OUTPUT_ROOT/status"

common_env=(
  "PM_PUB_ROOT=$PROSTATE_DATA_ROOT"
  "PM_REP_BASE=$PROSTATE_DATA_ROOT/prostate_replicates"
  "PM_TRUTH_DIR=$PROSTATE_DATA_ROOT/in_silico_true"
  "PM_SVCFIT_OUT=$parts"
  "PM_REPS=$rep"
  "SVCFIT_RLIB=$SVCFIT_RLIB"
)

env "${common_env[@]}" \
  PM_CONDITIONS=3m19,3m28,3m37,3m46,3m55,3m64,3m73,3m82,3m91 \
  PM_OUT_NAME="svcfit_chrx_rep${rep}.tsv" \
  "$SVCFIT_RSCRIPT" "$PROSTATE_SCRIPT_DIR/06a_run_svcfit_chrx.R"

env "${common_env[@]}" \
  PM_CONDITIONS=4m,5m \
  PM_OUT_NAME="svcfit_45_rep${rep}.tsv" \
  "$SVCFIT_RSCRIPT" "$PROSTATE_SCRIPT_DIR/06a_run_svcfit_chrx.R"

for f in "$parts/svcfit_chrx_rep${rep}.tsv" "$parts/svcfit_45_rep${rep}.tsv"; do
  [[ -s "$f" ]] || { echo "ERROR: missing output $f" >&2; exit 1; }
done
printf 'replicate\t%s\nworkflow_commit\t%s\nsvcfit_commit\t%s\nstatus\tok\n' \
  "$rep" "$EXPECTED_WORKFLOW_COMMIT" "$EXPECTED_SVCFIT_COMMIT" \
  > "$PROSTATE_OUTPUT_ROOT/status/rep${rep}.tsv"
echo "Done: corrected prostate SVCFit replicate $rep"
