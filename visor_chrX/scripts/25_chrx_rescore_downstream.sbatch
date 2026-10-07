#!/bin/bash
#SBATCH --job-name=chrx_rescore_ds
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#
# After the rescore array: copy the unchanged SVclone score table beside the new SVCFit
# scores, then rerun the replicate bootstrap (19) and the SVCFit-vs-SVclone comparison (23)
# with CHRX_DIR pointed at the rescore root. Defaults (B=10000, BOOT_SEED=1) are unchanged.

set -euo pipefail
: "${CHRX_DIR:?}" "${CHRX_SCORE_ROOT:?}" "${CHRX_SCRIPTS:?}" "${COV:?}" "${RSCRIPT_BIN:?}"
: "${WORKFLOW_REPO:?}" "${EXPECTED_WORKFLOW_COMMIT:?}" "${N_REPS:?}"

[[ "$(git -C "$WORKFLOW_REPO" rev-parse HEAD)" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
  echo "ERROR: workflow checkout moved after submission" >&2; exit 1; }

for r in $(seq 1 "$N_REPS"); do
  dst="${CHRX_SCORE_ROOT}/scoring_rep${r}"
  [[ -s "$dst/chrx_svcfit_scores_c${COV}.tsv" ]] || { echo "ERROR: missing rescore for rep ${r}" >&2; exit 1; }
  for f in chrx_svclone_scores_c${COV}.tsv; do
    src="${CHRX_DIR}/scoring_rep${r}/$f"
    [[ -s "$src" ]] || { echo "ERROR: missing $src" >&2; exit 1; }
    cp -p "$src" "$dst/$f"
  done
done
mkdir -p "${CHRX_SCORE_ROOT}/scoring"

export COV CHRX_DIR="$CHRX_SCORE_ROOT"
"$RSCRIPT_BIN" "${CHRX_SCRIPTS}/19_chrx_bootstrap.R"
"$RSCRIPT_BIN" "${CHRX_SCRIPTS}/23_chrx_compare_bootstrap.R"
echo "DONE $(date -u +%FT%TZ)"
