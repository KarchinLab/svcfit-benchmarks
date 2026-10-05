#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
manifest=""; output_root=""; rscript=""; package=""; truth=""; rlib=""
array="0,2,6,12%4"; submit=0; force=0
usage(){ echo "Usage: $0 --manifest FILE --output-root DIR --rscript FILE --svcfit-package DIR --truth-dir DIR --rlib DIR [--array SPEC] [--force] [--submit]"; }
while (($#)); do case "$1" in
  --manifest) manifest=${2:?};shift 2;; --output-root) output_root=${2:?};shift 2;;
  --rscript) rscript=${2:?};shift 2;; --svcfit-package) package=${2:?};shift 2;;
  --truth-dir) truth=${2:?};shift 2;; --rlib) rlib=${2:?};shift 2;;
  --array) array=${2:?};shift 2;; --force) force=1;shift;; --submit) submit=1;shift;;
  -h|--help) usage;exit 0;; *) usage >&2;exit 2;; esac; done
for p in "$manifest" "$output_root" "$rscript" "$package" "$truth" "$rlib"; do
  [[ "$p" = /* ]] || { echo "ERROR: absolute paths required" >&2;exit 2; }
done
[[ -s "$manifest" && -x "$rscript" && -d "$package" && -d "$truth" && -d "$rlib" ]] || exit 2
: "${EXPECTED_WORKFLOW_COMMIT:?Set EXPECTED_WORKFLOW_COMMIT to the audited workflow commit}"
: "${EXPECTED_SVCFIT_COMMIT:?Set EXPECTED_SVCFIT_COMMIT to the audited SVCFit commit}"
[[ "$EXPECTED_WORKFLOW_COMMIT" =~ ^[0-9a-f]{40}$ && "$EXPECTED_SVCFIT_COMMIT" =~ ^[0-9a-f]{40}$ ]] || {
  echo "ERROR: expected commits must contain 40 lowercase hexadecimal characters" >&2; exit 2;
}
workflow_repo=$(git -C "$script_dir" rev-parse --show-toplevel)
workflow_commit=$(git -C "$workflow_repo" rev-parse HEAD)
svcfit_commit=$(git -C "$package" rev-parse HEAD)
[[ "$workflow_commit" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
  echo "ERROR: workflow commit is $workflow_commit; expected $EXPECTED_WORKFLOW_COMMIT" >&2; exit 3;
}
[[ "$svcfit_commit" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
  echo "ERROR: SVCFit commit is $svcfit_commit; expected $EXPECTED_SVCFIT_COMMIT" >&2; exit 3;
}
if ((submit)); then
  [[ -z "$(git -C "$workflow_repo" status --porcelain)" ]] || { echo "ERROR: workflow checkout is dirty" >&2; exit 3; }
  [[ -z "$(git -C "$package" status --porcelain)" ]] || { echo "ERROR: SVCFit checkout is dirty" >&2; exit 3; }
fi
rows=$(awk 'END{print NR-1}' "$manifest"); ((rows==2250)) || { echo "ERROR: manifest rows=$rows" >&2;exit 3; }
echo "Preflight: manifest=2250 conditions; array=$array; output=$output_root"
echo "Workflow commit: $workflow_commit"
echo "SVCFit commit:   $svcfit_commit"
((submit)) || { echo "DRY RUN ONLY";exit 0; }
mkdir -p "$output_root/log"
job=$(sbatch --parsable --array="$array" \
  --output="$output_root/log/svcfit_3arm_%A_%a.out" --error="$output_root/log/svcfit_3arm_%A_%a.err" \
  --export="ALL,THREE_ARM_MANIFEST=$manifest,THREE_ARM_OUTPUT_ROOT=$output_root,VISOR_PIPELINE_DIR=$script_dir,WORKFLOW_REPO=$workflow_repo,EXPECTED_WORKFLOW_COMMIT=$workflow_commit,SVCFIT_RSCRIPT=$rscript,SVCFIT_PKG_DIR=$package,EXPECTED_SVCFIT_COMMIT=$svcfit_commit,SVCFIT_TRUTH_DIR=$truth,SVCFIT_RLIB=$rlib,FORCE_SVCFIT_EXTRACT=$force" \
  "$script_dir/08_svcfit_three_arm_array.sh")
printf 'job=%s\nworkflow_commit=%s\nsvcfit_commit=%s\narray=%s\nmanifest=%s\n' \
  "$job" "$workflow_commit" "$svcfit_commit" "$array" "$manifest" \
  > "$output_root/SVCFIT-SUBMISSION-${job}.txt"
echo "Submitted SVCFit three-arm extraction job $job"
