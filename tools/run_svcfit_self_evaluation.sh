#!/usr/bin/env bash
# Submit the canonical SVCFit software self-evaluation. Dry-run is the default.

set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage: run_svcfit_self_evaluation.sh [--dry-run | --submit] [--run-id ID]

Required environment:
  EXPECTED_WORKFLOW_COMMIT  Audited 40-character svcfit_workflows commit.

The SVCFit source checkout is pinned to commit
bedf5efc334841f854c90f5e489c86b854a547a2. The default is a no-write dry-run;
--submit is required to create a run directory and submit the Slurm job.
USAGE
}

submit=false
run_id=
while (($#)); do
    case "$1" in
        --dry-run) submit=false; shift ;;
        --submit) submit=true; shift ;;
        --run-id) run_id=${2-}; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) printf 'ERROR: unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

readonly script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly workflow_repo="$(git -C "$script_dir" rev-parse --show-toplevel)"
readonly config_file="${VISOR_CONFIG:-$workflow_repo/config.local.sh}"
readonly expected_svcfit_commit=bedf5efc334841f854c90f5e489c86b854a547a2

: "${EXPECTED_WORKFLOW_COMMIT:?Set EXPECTED_WORKFLOW_COMMIT to the audited full workflow commit}"
[[ "$EXPECTED_WORKFLOW_COMMIT" =~ ^[0-9a-f]{40}$ ]] || {
    printf 'ERROR: EXPECTED_WORKFLOW_COMMIT must contain 40 lowercase hex characters\n' >&2
    exit 1
}
[[ -r "$config_file" ]] || { printf 'ERROR: missing configuration: %s\n' "$config_file" >&2; exit 1; }
# shellcheck disable=SC1090
source "$config_file"

readonly svcfit_repo="${SVCFIT_PKG_DIR:?Set SVCFIT_PKG_DIR in the site configuration}"
readonly rscript_bin="${RSCRIPT:?Set RSCRIPT in the site configuration}"
readonly svcfit_env="${ENV_SVCFIT:?Set ENV_SVCFIT in the site configuration}"
readonly base_r_lib="${RLIB:?Set RLIB in the site configuration}"
readonly reticulate_python="${RETICULATE_PYTHON:?Set RETICULATE_PYTHON in the site configuration}"
readonly runtime_loader="${SVCFIT_RUNTIME_LOADER:?Set SVCFIT_RUNTIME_LOADER in the site configuration}"

if [[ -n "${SBATCH_BIN:-}" ]]; then
    sbatch_bin=$SBATCH_BIN
elif [[ -n "${SLURM_BIN:-}" ]]; then
    sbatch_bin=$SLURM_BIN/sbatch
else
    sbatch_bin=$(command -v sbatch || true)
fi
readonly sbatch_bin

for path in "$workflow_repo" "$svcfit_repo" "$svcfit_env"; do
    [[ -d "$path" ]] || { printf 'ERROR: required directory is absent: %s\n' "$path" >&2; exit 1; }
done
for path in "$rscript_bin" "$svcfit_env/bin/R" "$reticulate_python" "$sbatch_bin"; do
    [[ -x "$path" ]] || { printf 'ERROR: required executable is absent: %s\n' "$path" >&2; exit 1; }
done
[[ -r "$runtime_loader" ]] || { printf 'ERROR: runtime loader is absent: %s\n' "$runtime_loader" >&2; exit 1; }

workflow_commit=$(git -C "$workflow_repo" rev-parse HEAD)
svcfit_commit=$(git -C "$svcfit_repo" rev-parse HEAD)
[[ "$workflow_commit" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
    printf 'ERROR: workflow commit is %s; expected %s\n' "$workflow_commit" "$EXPECTED_WORKFLOW_COMMIT" >&2
    exit 1
}
[[ "$svcfit_commit" == "$expected_svcfit_commit" ]] || {
    printf 'ERROR: SVCFit commit is %s; expected %s\n' "$svcfit_commit" "$expected_svcfit_commit" >&2
    exit 1
}
if [[ -n "$(git -C "$workflow_repo" status --porcelain)" ]]; then
    if $submit; then
        printf 'ERROR: workflow checkout is dirty\n' >&2
        exit 1
    fi
    printf 'WARNING: workflow checkout is dirty; submission would be rejected\n' >&2
fi
[[ -z "$(git -C "$svcfit_repo" status --porcelain)" ]] || {
    printf 'ERROR: SVCFit checkout is dirty\n' >&2
    exit 1
}

run_id=${run_id:-package-tests-${svcfit_commit:0:7}-$(date -u '+%Y%m%dT%H%M%SZ')}
[[ "$run_id" =~ ^[A-Za-z0-9._-]+$ ]] || { printf 'ERROR: invalid --run-id: %s\n' "$run_id" >&2; exit 2; }
readonly run_root="${RUN_ROOT:-$ANALYSIS_ROOT/svcfit_self_evaluation/runs/$run_id}"

printf 'Evaluation scope: installed SVCFit package tests (software behavior)\n'
printf 'Run root:        %s\n' "$run_root"
printf 'Workflow commit: %s\n' "$workflow_commit"
printf 'SVCFit commit:   %s\n' "$svcfit_commit"
if ! $submit; then
    printf 'DRY RUN: no directory created and no job submitted.\n'
    exit 0
fi

[[ ! -e "$run_root" ]] || { printf 'ERROR: RUN_ROOT already exists: %s\n' "$run_root" >&2; exit 1; }
mkdir -p "$run_root/log" "$run_root/provenance" "$run_root/library" \
    "$run_root/cache" "$run_root/r-user" "$run_root/tmp"
printf '%s\n' "$workflow_commit" > "$run_root/provenance/svcfit_workflows.commit.txt"
printf '%s\n' "$svcfit_commit" > "$run_root/provenance/SVCFit.commit.txt"
sha256sum "$config_file" > "$run_root/provenance/site-config.sha256"
cp -p "$runtime_loader" "$run_root/provenance/site-runtime.sh"
if [[ -r "${SVCFIT_RUNTIME_LOCKFILE:-}" ]]; then
    cp -p "$SVCFIT_RUNTIME_LOCKFILE" "$run_root/provenance/conda-explicit-linux-64.txt"
fi

slurm_site_args=()
[[ -n "${SLURM_ACCOUNT:-}" ]] && slurm_site_args+=(--account="$SLURM_ACCOUNT")
[[ -n "${SLURM_PARTITION:-}" ]] && slurm_site_args+=(--partition="$SLURM_PARTITION")
[[ -n "${SLURM_QOS:-}" ]] && slurm_site_args+=(--qos="$SLURM_QOS")

export RUN_ROOT="$run_root" RUN_LIBRARY="$run_root/library" WORKFLOW_REPO="$workflow_repo"
export EXPECTED_WORKFLOW_COMMIT="$workflow_commit"
export EXPECTED_SVCFIT_COMMIT="$svcfit_commit" SVCFIT_PKG_DIR="$svcfit_repo"
export RSCRIPT="$rscript_bin" ENV_SVCFIT="$svcfit_env" RLIB="$base_r_lib"
export RETICULATE_PYTHON="$reticulate_python"

job_id=$("$sbatch_bin" --parsable "${slurm_site_args[@]}" \
    --chdir="$workflow_repo" --job-name=svcfit-self-eval \
    --mem="${SELF_EVAL_MEM:-8G}" --time="${SELF_EVAL_TIME:-00:30:00}" \
    --output="$run_root/log/self-evaluation-%j.out" \
    --error="$run_root/log/self-evaluation-%j.err" --export=ALL \
    "$script_dir/svcfit_self_evaluation.sbatch")
job_id=${job_id%%;*}

{
    printf 'job=%s\n' "$job_id"
    printf 'submitted=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf 'scope=installed-package-tests\n'
    printf 'workflow_commit=%s\n' "$workflow_commit"
    printf 'svcfit_commit=%s\n' "$svcfit_commit"
    printf 'run_root=%s\n' "$run_root"
} > "$run_root/SUBMISSION-${job_id}.txt"

printf 'Submitted SVCFit software self-evaluation: %s\n' "$job_id"
printf 'Submission record: %s\n' "$run_root/SUBMISSION-${job_id}.txt"
