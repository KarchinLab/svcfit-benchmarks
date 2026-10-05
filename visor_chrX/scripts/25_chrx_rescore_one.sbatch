#!/bin/bash
#SBATCH --job-name=chrx_rescore
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#
# Rescore ONE replicate (REP = array task id) with a pinned SVCFit build: step 13 only.
# Calls, segmentation and the 09 join are read from CHRX_DIR unchanged; scores go to
# CHRX_SCORE_ROOT. Submitted by 25_chrx_rescore_svcfit.sh, which sets every variable below.

set -euo pipefail
: "${CHRX_DIR:?}" "${CHRX_SCORE_ROOT:?}" "${CHRX_SCRIPTS:?}" "${COV:?}" "${RSCRIPT_BIN:?}"
: "${SVCFIT_R_LIB:?}" "${RLIB:?}" "${WORKFLOW_REPO:?}" "${EXPECTED_WORKFLOW_COMMIT:?}"
: "${EXPECTED_SVCFIT_COMMIT:?}" "${SVCFIT_R_LIB_COMMIT_FILE:?}" "${SLURM_ARRAY_TASK_ID:?}"

[[ "$(git -C "$WORKFLOW_REPO" rev-parse HEAD)" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
  echo "ERROR: workflow checkout moved after submission" >&2; exit 1; }
[[ "$(tr -d '[:space:]' < "$SVCFIT_R_LIB_COMMIT_FILE")" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
  echo "ERROR: SVCFit library provenance does not match $EXPECTED_SVCFIT_COMMIT" >&2; exit 1; }

export REP="$SLURM_ARRAY_TASK_ID" COV CHRX_DIR CHRX_SCORE_ROOT
export R_LIBS="${SVCFIT_R_LIB}:${RLIB}"
"$RSCRIPT_BIN" "${CHRX_SCRIPTS}/13_chrx_score_svcfit.R"

out="${CHRX_SCORE_ROOT}/scoring_rep${REP}/chrx_svcfit_scores_c${COV}.tsv"
[[ -s "$out" ]] || { echo "ERROR: replicate ${REP} produced no $out" >&2; exit 1; }
echo "RESULT: replicate ${REP} rescored -> $out ($(( $(wc -l < "$out") - 1 )) rows)"
