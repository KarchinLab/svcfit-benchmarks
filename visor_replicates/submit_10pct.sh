#!/bin/bash
# Submit VISOR bootstrap pipeline for p10 purity only.
# 450 tasks = 30 reps x 5 exps x 3 conditions (p10m10, p10m30, p10m50)
#
# Strategy: pipe each existing pipeline script through sed to swap
# pipeline_common.sh → pipeline_common_10pct.sh, then submit to SLURM.
# --array=0-449%10 on the command line overrides the #SBATCH directive in each script.
# This does NOT touch or interfere with already-running jobs.
#
# Dependency chain mirrors submit_all.sh:
#   00 → 01 → 02 → 03  (04 is not submitted)

set -euo pipefail

# --- machine-specific paths ---------------------------------------------------
# Absolute paths below come from config.local.sh, which is per-machine and never
# committed. Located by walking UP from this script. See README.md.
if [[ -z "${VISOR_CONFIG_LOADED:-}" ]]; then
    _cfg="${VISOR_CONFIG:-}"
    if [[ -z "$_cfg" ]]; then
        _d="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        while [[ "$_d" != "/" ]]; do
            [[ -r "$_d/config.local.sh" ]] && { _cfg="$_d/config.local.sh"; break; }
            _d="$(dirname "$_d")"
        done
    fi
    [[ -n "$_cfg" && -r "$_cfg" ]] || {
        echo "ERROR: config.local.sh not found by walking up from ${BASH_SOURCE[0]}." >&2
        echo "       cp config.example.sh config.local.sh at the repo root and edit it." >&2
        echo "       See README.md." >&2
        exit 1
    }
    # shellcheck disable=SC1090
    source "$_cfg"
    unset _cfg _d
fi

cd "$(dirname "$0")"

mkdir -p log

SCRIPT_DIR="$(pwd)"

# Submit a pipeline script with pipeline_common_10pct.sh swapped in.
# Additional sbatch args (e.g. --dependency) are passed as extra arguments.
submit_10pct() {
    local script="${SCRIPT_DIR}/$1"
    shift
    sed 's|pipeline_common\.sh|pipeline_common_10pct.sh|g' "$script" \
        | sbatch --parsable --array=0-449%10 "$@"
}

echo "Submitting VISOR p10 purity pipeline (450 tasks each)..."

JOB0=$(submit_10pct 00_visor_shorts.sh)
echo "00_visor_shorts    : $JOB0  (array 0-449)"

JOB1=$(submit_10pct 01_manta_svtyper.sh --dependency=afterok:${JOB0})
echo "01_manta_svtyper   : $JOB1  (array 0-449, after $JOB0)"

JOB2=$(submit_10pct 02_snp_pipeline.sh --dependency=aftercorr:${JOB1})
echo "02_snp_pipeline    : $JOB2  (array 0-449, after each task of $JOB1)"

JOB3=$(submit_10pct 03_facet.sh --dependency=afterok:${JOB2})
echo "03_facet           : $JOB3  (array 0-449, after $JOB2)"

echo ""
echo "Monitor with: squeue -u \$USER"
echo "After job 03 completes, run the SVclone arms and the SVCFit three-arm stages (see README.md)."
