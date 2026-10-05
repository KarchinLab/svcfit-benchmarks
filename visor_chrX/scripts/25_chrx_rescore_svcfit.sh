#!/usr/bin/env bash
# Rescore the chrX simulation with a pinned, installed SVCFit build, without rerunning the
# simulation, calling, segmentation or the 09 join, and without touching the accepted
# scoring_rep* tables. Dry-run by default.
#
# Required environment:
#   EXPECTED_WORKFLOW_COMMIT   40-character svcfit_workflows commit
#   EXPECTED_SVCFIT_COMMIT     40-character SVCFit commit the library was built from
#   SVCFIT_R_LIB               R library holding that SVCFit build
#   SVCFIT_R_LIB_COMMIT_FILE   file naming the commit SVCFIT_R_LIB was built from
# Optional: RUN_ROOT, N_REPS (default 30), COV (default 50).
#
# Usage: 25_chrx_rescore_svcfit.sh [--submit]

set -Eeuo pipefail
submit=false
case "${1:-}" in --submit) submit=true ;; ""|--dry-run) ;; *) echo "Usage: $0 [--submit]" >&2; exit 2 ;; esac

readonly script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly workflow_repo="$(git -C "$script_dir" rev-parse --show-toplevel)"
# shellcheck disable=SC1091
source "${VISOR_CONFIG:-$workflow_repo/config.local.sh}"
: "${CHRX_DIR:?Set CHRX_DIR in the site configuration}" "${RLIB:?}" "${SVCFIT_R:?}"
: "${EXPECTED_WORKFLOW_COMMIT:?}" "${EXPECTED_SVCFIT_COMMIT:?}"
: "${SVCFIT_R_LIB:?}" "${SVCFIT_R_LIB_COMMIT_FILE:?}"
readonly cov=${COV:-50} n_reps=${N_REPS:-30}
sbatch_bin="${SBATCH_BIN:-${SLURM_BIN:+$SLURM_BIN/sbatch}}"; sbatch_bin="${sbatch_bin:-$(command -v sbatch || true)}"

workflow_commit=$(git -C "$workflow_repo" rev-parse HEAD)
[[ "$workflow_commit" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
  echo "ERROR: workflow commit is $workflow_commit; expected $EXPECTED_WORKFLOW_COMMIT" >&2; exit 1; }
lib_commit=$(tr -d '[:space:]' < "$SVCFIT_R_LIB_COMMIT_FILE")
[[ "$lib_commit" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
  echo "ERROR: SVCFIT_R_LIB was built from $lib_commit; expected $EXPECTED_SVCFIT_COMMIT" >&2; exit 1; }
loaded=$(R_LIBS="$SVCFIT_R_LIB:$RLIB" "$SVCFIT_R" --vanilla -e \
  'suppressMessages(library(SVCFit)); h <- c("hemizygous_del_svcf","hemizygous_dup_svcf","resolve_hemizygous_svcf"); stopifnot(all(h %in% ls("package:SVCFit"))); cat(find.package("SVCFit"))' 2>/dev/null) || {
  echo "ERROR: R cannot load SVCFit with the hemizygous helpers from $SVCFIT_R_LIB" >&2; exit 1; }
[[ "$loaded" == "$(cd "$SVCFIT_R_LIB" && pwd -P)/SVCFit" ]] || {
  echo "ERROR: R would load SVCFit from $loaded, not SVCFIT_R_LIB" >&2; exit 1; }
for r in $(seq 1 "$n_reps"); do
  for d in "calls_rep$r" "cn_bar10k_m_rep$r" "scoring_rep$r"; do
    [[ -d "$CHRX_DIR/$d" ]] || { echo "ERROR: missing input $CHRX_DIR/$d" >&2; exit 1; }
  done
done
[[ -s "$CHRX_DIR/truth/chrx_sv_ground_truth_c${cov}.tsv" ]] || { echo "ERROR: missing truth table" >&2; exit 1; }

readonly run_root=${RUN_ROOT:-$(dirname "$CHRX_DIR")/rescore-svcfit-${EXPECTED_SVCFIT_COMMIT:0:7}-$(date -u +%Y%m%dT%H%MZ)}
[[ ! -e "$run_root" ]] || { echo "ERROR: RUN_ROOT already exists: $run_root" >&2; exit 1; }

printf 'Input (read-only): %s\n' "$CHRX_DIR"
printf 'Run root:          %s\n' "$run_root"
printf 'Replicates:        1-%s (COV=%s)\n' "$n_reps" "$cov"
printf 'Workflow commit:   %s\n' "$workflow_commit"
printf 'SVCFit library:    %s (%s)\n' "$loaded" "$lib_commit"
$submit || { echo "DRY RUN: no directory created and no job submitted."; exit 0; }
[[ -z "$(git -C "$workflow_repo" status --porcelain)" ]] || { echo "ERROR: workflow checkout is dirty" >&2; exit 1; }

mkdir -p "$run_root/log" "$run_root/provenance"
printf '%s\n' "$workflow_commit" > "$run_root/provenance/svcfit_workflows.commit.txt"
printf '%s\n' "$lib_commit" > "$run_root/provenance/SVCFit.commit.txt"
printf 'r_libs\tloaded_from\tinstalled_commit\n%s\t%s\t%s\n' "$SVCFIT_R_LIB:$RLIB" "$loaded" "$lib_commit" \
  > "$run_root/provenance/svcfit-r-library.tsv"
export CHRX_DIR CHRX_SCORE_ROOT="$run_root" CHRX_SCRIPTS="$script_dir" COV="$cov" N_REPS="$n_reps" \
  RSCRIPT_BIN="$SVCFIT_R" SVCFIT_R_LIB RLIB WORKFLOW_REPO="$workflow_repo" \
  EXPECTED_WORKFLOW_COMMIT EXPECTED_SVCFIT_COMMIT SVCFIT_R_LIB_COMMIT_FILE
site=(); [[ -n "${SLURM_PARTITION:-}" ]] && site+=(--partition="$SLURM_PARTITION")
score_job=$("$sbatch_bin" --parsable "${site[@]}" --array="1-${n_reps}%15" --export=ALL \
  --output="$run_root/log/rescore_%A_%a.out" --error="$run_root/log/rescore_%A_%a.err" \
  "$script_dir/25_chrx_rescore_one.sbatch")
ds_job=$("$sbatch_bin" --parsable "${site[@]}" --dependency="afterok:${score_job}" --export=ALL \
  --output="$run_root/log/downstream_%j.out" --error="$run_root/log/downstream_%j.err" \
  "$script_dir/25_chrx_rescore_downstream.sbatch")
printf 'score_job\t%s\ndownstream_job\t%s\n' "$score_job" "$ds_job" > "$run_root/provenance/submitted-jobs.tsv"
echo "Submitted rescore array $score_job and downstream job $ds_job"
