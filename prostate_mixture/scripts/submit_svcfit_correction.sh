#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
output_root=""; data_root=""; rscript=""; package=""; rlib=""; rlib_commit_file=""
array="0-29%10"; submit=0
usage() {
  echo "Usage: $0 --output-root DIR --data-root DIR --rscript FILE --svcfit-package DIR --rlib DIR --rlib-commit-file FILE [--array SPEC] [--submit]"
}
while (($#)); do
  case "$1" in
    --output-root) output_root=${2:?}; shift 2;;
    --data-root) data_root=${2:?}; shift 2;;
    --rscript) rscript=${2:?}; shift 2;;
    --svcfit-package) package=${2:?}; shift 2;;
    --rlib) rlib=${2:?}; shift 2;;
    --rlib-commit-file) rlib_commit_file=${2:?}; shift 2;;
    --array) array=${2:?}; shift 2;;
    --submit) submit=1; shift;;
    -h|--help) usage; exit 0;;
    *) usage >&2; exit 2;;
  esac
done
for p in "$output_root" "$data_root" "$rscript" "$package" "$rlib" "$rlib_commit_file"; do
  [[ "$p" = /* ]] || { echo "ERROR: absolute paths required" >&2; exit 2; }
done
[[ -d "$data_root/prostate_replicates" && -d "$data_root/in_silico_true" ]] || exit 2
[[ -x "$rscript" && -d "$package" && -d "$rlib/SVCFit" && -s "$rlib_commit_file" ]] || exit 2

: "${EXPECTED_WORKFLOW_COMMIT:?Set EXPECTED_WORKFLOW_COMMIT to the audited workflow commit}"
: "${EXPECTED_SVCFIT_COMMIT:?Set EXPECTED_SVCFIT_COMMIT to the audited SVCFit commit}"
[[ "$EXPECTED_WORKFLOW_COMMIT" =~ ^[0-9a-f]{40}$ && "$EXPECTED_SVCFIT_COMMIT" =~ ^[0-9a-f]{40}$ ]] || exit 2
workflow_repo=$(git -C "$script_dir" rev-parse --show-toplevel)
workflow_commit=$(git -C "$workflow_repo" rev-parse HEAD)
svcfit_commit=$(git -C "$package" rev-parse HEAD)
rlib_commit=$(tr -d '[:space:]' < "$rlib_commit_file")
[[ "$workflow_commit" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
  echo "ERROR: workflow commit is $workflow_commit; expected $EXPECTED_WORKFLOW_COMMIT" >&2; exit 3;
}
[[ "$svcfit_commit" == "$EXPECTED_SVCFIT_COMMIT" && "$rlib_commit" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
  echo "ERROR: SVCFit source/library provenance does not match $EXPECTED_SVCFIT_COMMIT" >&2; exit 3;
}
if ((submit)); then
  [[ -z "$(git -C "$workflow_repo" status --porcelain)" ]] || { echo "ERROR: workflow checkout is dirty" >&2; exit 3; }
  [[ -z "$(git -C "$package" status --porcelain)" ]] || { echo "ERROR: SVCFit checkout is dirty" >&2; exit 3; }
fi

echo "Preflight: prostate replicates=30; array=$array; output=$output_root"
echo "Workflow commit: $workflow_commit"
echo "SVCFit commit:   $svcfit_commit"
((submit)) || { echo "DRY RUN ONLY"; exit 0; }
mkdir -p "$output_root/log"
job=$(sbatch --parsable --array="$array" \
  --output="$output_root/log/prostate_svcf_%A_%a.out" \
  --error="$output_root/log/prostate_svcf_%A_%a.err" \
  --export="ALL,PROSTATE_OUTPUT_ROOT=$output_root,PROSTATE_DATA_ROOT=$data_root,PROSTATE_SCRIPT_DIR=$script_dir,WORKFLOW_REPO=$workflow_repo,EXPECTED_WORKFLOW_COMMIT=$workflow_commit,SVCFIT_PKG_DIR=$package,EXPECTED_SVCFIT_COMMIT=$svcfit_commit,SVCFIT_RSCRIPT=$rscript,SVCFIT_RLIB=$rlib,SVCFIT_RLIB_COMMIT_FILE=$rlib_commit_file" \
  "$script_dir/09_svcfit_correction_array.sh")
printf 'job=%s\nworkflow_commit=%s\nsvcfit_commit=%s\narray=%s\ndata_root=%s\nrlib=%s\n' \
  "$job" "$workflow_commit" "$svcfit_commit" "$array" "$data_root" "$rlib" \
  > "$output_root/SUBMISSION-${job}.txt"
echo "Submitted corrected prostate SVCFit job $job"
