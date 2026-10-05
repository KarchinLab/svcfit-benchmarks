#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source_root=""; result_root=""; purity_table=""; conda_sh=""; env_dir=""
array="0,3,6,12,14%5"; submit=0; force=0
usage(){ echo "Usage: $0 --source-root DIR --result-root DIR --purity-table CSV --svclone-conda-sh FILE --svclone-env DIR [--array SPEC] [--force] [--submit]"; }
while (($#)); do case "$1" in
  --source-root) source_root=${2:?};shift 2;; --result-root) result_root=${2:?};shift 2;;
  --purity-table) purity_table=${2:?};shift 2;; --svclone-conda-sh) conda_sh=${2:?};shift 2;;
  --svclone-env) env_dir=${2:?};shift 2;; --array) array=${2:?};shift 2;;
  --force) force=1;shift;; --submit) submit=1;shift;; -h|--help) usage;exit 0;; *) usage >&2;exit 2;; esac; done
for p in "$source_root" "$result_root" "$purity_table" "$conda_sh" "$env_dir"; do [[ "$p" = /* ]] || { echo "ERROR: absolute paths required" >&2;exit 2;}; done
[[ -d "$source_root" && -r "$purity_table" && -r "$conda_sh" && -d "$env_dir" && "$source_root" != "$result_root" ]] || exit 2

rows=$(awk 'END{print NR-1}' "$purity_table"); usable=$(awk -F, 'NR>1&&$4!="NA"{n++}END{print n+0}' "$purity_table")
((rows==2250 && usable==1727)) || { echo "ERROR: purity table rows=$rows usable=$usable" >&2;exit 3;}
echo "Preflight: 2250 conditions; 1727 FACETS-usable; 523 expected failures"
echo "Array: $array; calling disabled; result root: $result_root"
((submit)) || { echo "DRY RUN ONLY";exit 0;}
mkdir -p "$result_root/log"
job=$(sbatch --parsable --array="$array" --output="$result_root/log/fair_%A_%a.out" --error="$result_root/log/fair_%A_%a.err" \
  --export="ALL,SOURCE_REPLICATES_DIR=$source_root,FAIR_RESULT_ROOT=$result_root,FACETS_PURITY_TABLE=$purity_table,VISOR_PIPELINE_DIR=$script_dir,VISOR_SVCLONE_CONDA_SH=$conda_sh,VISOR_SVCLONE_ENV=$env_dir,FORCE_FAIR_RERUN=$force" \
  "$script_dir/05_svclone_fair_from_existing_calls.sh")
printf 'job=%s\ncommit=%s\narray=%s\ncalling=disabled\n' "$job" "$(git -C "$script_dir" rev-parse HEAD)" "$array" > "$result_root/SUBMISSION-${job}.txt"
echo "Submitted fair-arm job $job"
