#!/usr/bin/env bash

# Submit SVclone for all 30 replicates x 5 experiments x 15 conditions using
# existing SV calls. Dry-run is the default; --submit is intentionally required.

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source_root=""
result_root=""
config=""
svclone_conda_sh=""
svclone_env=""
array_spec="0-2249%25"
submit=0
force=0

usage() {
    echo "Usage: $0 --source-root DIR --result-root DIR --config FILE --svclone-conda-sh FILE --svclone-env DIR [--array SPEC] [--force] [--submit]"
}

while (( $# )); do
    case "$1" in
        --source-root) source_root=${2:?}; shift 2 ;;
        --result-root) result_root=${2:?}; shift 2 ;;
        --config) config=${2:?}; shift 2 ;;
        --svclone-conda-sh) svclone_conda_sh=${2:?}; shift 2 ;;
        --svclone-env) svclone_env=${2:?}; shift 2 ;;
        --array) array_spec=${2:?}; shift 2 ;;
        --force) force=1; shift ;;
        --submit) submit=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

source_root=${source_root%/}
result_root=${result_root%/}
[[ "$source_root" = /* && -d "$source_root" ]] || {
    echo "ERROR: --source-root must be an existing absolute directory" >&2; exit 2;
}
[[ "$result_root" = /* && "$result_root" != "$source_root" ]] || {
    echo "ERROR: --result-root must be absolute and differ from --source-root" >&2; exit 2;
}
[[ "$config" = /* && -r "$config" ]] || {
    echo "ERROR: --config must be an absolute readable file" >&2; exit 2;
}
[[ "$svclone_conda_sh" = /* && -r "$svclone_conda_sh" ]] || {
    echo "ERROR: --svclone-conda-sh must be an absolute readable file" >&2; exit 2;
}
[[ "$svclone_env" = /* && -d "$svclone_env" ]] || {
    echo "ERROR: --svclone-env must be an existing absolute directory" >&2; exit 2;
}
[[ "$array_spec" =~ ^[0-9]+-[0-9]+%[0-9]+$ ]] || {
    echo "ERROR: --array must look like 0-2249%25" >&2; exit 2;
}

purities=(10 10 10 20 20 20 40 40 40 60 60 60 80 80 80)
mixtures=(10 30 50 10 30 50 10 30 50 10 30 50 10 30 50)
missing=0
for task_id in $(seq 0 2249); do
    rep_idx=$(( task_id / 75 ))
    remainder=$(( task_id % 75 ))
    exp_idx=$(( remainder / 15 ))
    cond_idx=$(( remainder % 15 ))
    condition=c50p${purities[$cond_idx]}m${mixtures[$cond_idx]}
    work="$source_root/rep$((rep_idx + 1))/exp$((exp_idx + 1))/$condition"
    for path in \
        "$work/manta/results/variants/${condition}.vcf" \
        "$work/svtyp/svt_${condition}.vcf" \
        "$work/short/short.out/sim.srt.bam" \
        "$work/facet/${condition}.RData" \
        "$work/facet/${condition}.bed"; do
        if [[ ! -s "$path" ]]; then
            echo "MISSING: $path" >&2
            missing=$((missing + 1))
        fi
    done
done
(( missing == 0 )) || { echo "ERROR: preflight found $missing missing/empty inputs" >&2; exit 3; }

echo "Preflight passed: 2250 samples (30 replicates x 5 experiments x 15 conditions)"
echo "Source calls : $source_root"
echo "Result root  : $result_root"
echo "Array        : $array_spec"
echo "Calling      : disabled (existing VCFs only)"
echo "SVclone env  : $svclone_env"

if (( submit == 0 )); then
    echo "DRY RUN ONLY: add --submit to create the result root and submit SVclone"
    exit 0
fi

[[ ! -e "$result_root" || -d "$result_root" ]] || {
    echo "ERROR: result root exists and is not a directory: $result_root" >&2; exit 4;
}
mkdir -p "$result_root/log"

input_manifest="$result_root/INPUTS.tsv"
"$script_dir/build_input_manifest.sh" \
    --source-root "$source_root" \
    --result-root "$result_root" \
    --output "$input_manifest"
input_manifest_sha=$(sha256sum "$input_manifest" | awk '{print $1}')

job_id=$(sbatch --parsable \
    --array="$array_spec" \
    --output="$result_root/log/04_svclone_%A_%a.out" \
    --error="$result_root/log/04_svclone_%A_%a.err" \
    --export="ALL,SOURCE_REPLICATES_DIR=$source_root,VISOR_REPLICATES_DIR=$result_root,VISOR_INPUT_MANIFEST=$input_manifest,VISOR_CONFIG=$config,VISOR_PIPELINE_DIR=$script_dir,VISOR_SVCLONE_CONDA_SH=$svclone_conda_sh,VISOR_SVCLONE_ENV=$svclone_env,FORCE_DOWNSTREAM_RERUN=$force" \
    "$script_dir/04_svclone_from_existing_calls.sh")

submission_record="$result_root/SUBMISSION-${job_id}.txt"
printf '%s\n' \
    "Commit       : $(git -C "$script_dir" rev-parse HEAD)" \
    "Submitted    : $(date --iso-8601=seconds)" \
    "Job ID       : $job_id" \
    "Array        : $array_spec" \
    "Source calls : $source_root" \
    "Result root  : $result_root" \
    "Input manifest: $input_manifest" \
    "Manifest SHA : $input_manifest_sha" \
    "Calling      : disabled" \
    "SVclone env  : $svclone_env" \
    > "$submission_record"
cat "$submission_record" >> "$result_root/SUBMISSIONS.log"
printf '\n' >> "$result_root/SUBMISSIONS.log"

echo "Submitted downstream-only SVclone array: $job_id"
