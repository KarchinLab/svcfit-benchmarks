#!/usr/bin/env bash
# Submit the pinned S1 phylogeny benchmark (SVCFit e0d7e0b).

set -Eeuo pipefail

usage() {
  printf '%s\n' \
    'Usage: run_svcfit_and_evaluate.sh --scope smoke|full [--dry-run]' \
    '' \
    'Required environment: EXPECTED_WORKFLOW_COMMIT (full 40-character commit).' \
    'Optional: SVCFIT_R_LIB, an R library holding SVCFit installed from the' \
    'pinned commit, placed ahead of RLIB; SVCFIT_R_LIB_COMMIT_FILE must then name that commit.'
}

mode=e0d7e0b
scope=
dry_run=false
while (($#)); do
  case "$1" in
    --scope) scope=${2-}; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) printf 'ERROR: unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ "$scope" == smoke || "$scope" == full ]] || { printf 'ERROR: --scope is required\n' >&2; exit 2; }

readonly script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly workflow_repo="$(git -C "$script_dir" rev-parse --show-toplevel)"
readonly config_file="${VISOR_CONFIG:-$workflow_repo/config.local.sh}"
[[ -r "$config_file" ]] || { printf 'ERROR: missing configuration: %s\n' "$config_file" >&2; exit 1; }
# shellcheck disable=SC1090
source "$config_file"
readonly project_root="$PROJECT_ROOT"
readonly input_dir="${TREE_EVAL_LONGITUDINAL:-$project_root/02_data/tree_eval/longitudinal}"
readonly truth_dir="${TREE_EVAL_TRUTH_DIR:-$project_root/02_data/tree_eval/input_data/hack}"
readonly runtime_loader="${SVCFIT_RUNTIME_LOADER:?Set SVCFIT_RUNTIME_LOADER in the site configuration}"
if [[ -n "${SBATCH_BIN:-}" ]]; then
  sbatch_bin="$SBATCH_BIN"
elif [[ -n "${SLURM_BIN:-}" ]]; then
  sbatch_bin="$SLURM_BIN/sbatch"
else
  sbatch_bin="$(command -v sbatch || true)"
fi
readonly sbatch_bin
readonly corrected_commit=e0d7e0b8d704caa3cfb227bc3dbf1ee99d66fac7  # SVCFit commit used for the paper

: "${EXPECTED_WORKFLOW_COMMIT:?Set EXPECTED_WORKFLOW_COMMIT to the audited full workflow commit}"
[[ "$EXPECTED_WORKFLOW_COMMIT" =~ ^[0-9a-f]{40}$ ]] || { printf 'ERROR: EXPECTED_WORKFLOW_COMMIT must contain 40 lowercase hex characters\n' >&2; exit 1; }

readonly svcfit_repo="${SVCFIT_PKG_DIR:-$project_root/01_software/SVCFit}"
readonly expected_svcfit_commit=$corrected_commit
expected_no_tree=
expected_correct=
svcfit_r_lib="${RLIB:?Set RLIB in the site configuration}"
# Workers run library(SVCFit), not the checkout, so the installed copy must match the pin.
if [[ -n "${SVCFIT_R_LIB:-}" ]]; then
  [[ -d "$SVCFIT_R_LIB/SVCFit" ]] || { printf 'ERROR: SVCFIT_R_LIB has no SVCFit package: %s\n' "$SVCFIT_R_LIB" >&2; exit 1; }
  [[ -s "${SVCFIT_R_LIB_COMMIT_FILE:-}" ]] || { printf 'ERROR: SVCFIT_R_LIB_COMMIT_FILE is required with SVCFIT_R_LIB\n' >&2; exit 1; }
  svcfit_r_lib="$SVCFIT_R_LIB:$svcfit_r_lib"
fi

if [[ "$scope" == smoke ]]; then
  purities=(20)
  boots=(0)
  experiments=(exp1)
  array_spec=0
  expected_cases=1
else
  purities=(20 40 60 80)
  boots=(0 1 2 3 4)
  experiments=(exp1 exp2 exp3 exp4 exp5)
  array_spec=0-4
  expected_cases=100
fi

readonly run_id=${RUN_ID:-${mode}-${scope}-$(date -u '+%Y%m%dT%H%M%SZ')}
readonly run_root=${RUN_ROOT:-$project_root/03_analysis/tree_eval/runs/$run_id}
readonly output_dir="$run_root/output"
readonly eval_dir="$run_root/evaluation"
readonly log_dir="$run_root/log"
readonly provenance_dir="$run_root/provenance"

for path in "$workflow_repo" "$input_dir" "$truth_dir" "$svcfit_repo"; do
  [[ -d "$path" ]] || { printf 'ERROR: required directory is absent: %s\n' "$path" >&2; exit 1; }
done
[[ -r "$runtime_loader" ]] || { printf 'ERROR: runtime loader is absent: %s\n' "$runtime_loader" >&2; exit 1; }
[[ -x "$sbatch_bin" ]] || { printf 'ERROR: sbatch is absent: %s\n' "$sbatch_bin" >&2; exit 1; }
[[ ! -e "$run_root" ]] || { printf 'ERROR: RUN_ROOT already exists: %s\n' "$run_root" >&2; exit 1; }

workflow_commit=$(git -C "$workflow_repo" rev-parse HEAD)
svcfit_commit=$(git -C "$svcfit_repo" rev-parse HEAD)
[[ "$workflow_commit" == "$EXPECTED_WORKFLOW_COMMIT" ]] || { printf 'ERROR: workflow commit is %s; expected %s\n' "$workflow_commit" "$EXPECTED_WORKFLOW_COMMIT" >&2; exit 1; }
[[ "$svcfit_commit" == "$expected_svcfit_commit" ]] || { printf 'ERROR: SVCFit commit is %s; expected %s\n' "$svcfit_commit" "$expected_svcfit_commit" >&2; exit 1; }
if [[ -n "$(git -C "$workflow_repo" status --porcelain)" ]]; then
  if $dry_run; then
    printf 'WARNING: workflow checkout is dirty; submission would be rejected\n' >&2
  else
    printf 'ERROR: workflow checkout is dirty\n' >&2
    exit 1
  fi
fi
[[ -z "$(git -C "$svcfit_repo" status --porcelain)" ]] || { printf 'ERROR: SVCFit checkout is dirty\n' >&2; exit 1; }

svcfit_r_lib_commit=unverified
if [[ -n "${SVCFIT_R_LIB:-}" ]]; then
  svcfit_r_lib_commit=$(tr -d '[:space:]' < "$SVCFIT_R_LIB_COMMIT_FILE")
  [[ "$svcfit_r_lib_commit" == "$expected_svcfit_commit" ]] || { printf 'ERROR: SVCFIT_R_LIB was built from %s; expected %s\n' "$svcfit_r_lib_commit" "$expected_svcfit_commit" >&2; exit 1; }
fi
svcfit_loaded_from=$(R_LIBS_USER="$svcfit_r_lib" "${SVCFIT_R:?Set SVCFIT_R in the site configuration}" --vanilla -e 'cat(find.package("SVCFit"))' 2>/dev/null) || { printf 'ERROR: workers cannot find SVCFit on R_LIBS_USER=%s\n' "$svcfit_r_lib" >&2; exit 1; }
if [[ -n "${SVCFIT_R_LIB:-}" && "$svcfit_loaded_from" != "$(cd "$SVCFIT_R_LIB" && pwd -P)/SVCFit" ]]; then
  printf 'ERROR: workers would load SVCFit from %s, not SVCFIT_R_LIB\n' "$svcfit_loaded_from" >&2
  exit 1
fi
[[ "$svcfit_r_lib_commit" != unverified ]] || printf 'WARNING: installed SVCFit at %s is not verified against the pinned commit\n' "$svcfit_loaded_from" >&2

printf 'Mode:            %s\n' "$mode"
printf 'Scope:           %s\n' "$scope"
printf 'Run root:        %s\n' "$run_root"
printf 'Workflow commit: %s\n' "$workflow_commit"
printf 'SVCFit commit:   %s\n' "$svcfit_commit"
printf 'SVCFit library:  %s (%s)\n' "$svcfit_loaded_from" "$svcfit_r_lib_commit"
printf 'Expected cases:  %s\n' "$expected_cases"
if $dry_run; then
  printf 'DRY RUN: no directory created and no job submitted.\n'
  exit 0
fi

mkdir -p "$output_dir" "$eval_dir" "$log_dir" "$provenance_dir" "$run_root/runtime-home" "$run_root/cache" "$run_root/r-user"
printf '%s\n' "$workflow_commit" > "$provenance_dir/svcfit_workflows.commit.txt"
printf '%s\n' "$svcfit_commit" > "$provenance_dir/SVCFit.commit.txt"
{
  printf 'r_libs_user\tloaded_from\tinstalled_commit\n'
  printf '%s\t%s\t%s\n' "$svcfit_r_lib" "$svcfit_loaded_from" "$svcfit_r_lib_commit"
} > "$provenance_dir/svcfit-r-library.tsv"
cp -p "$runtime_loader" "$provenance_dir/site-runtime.sh"
if [[ -r "${SVCFIT_RUNTIME_LOCKFILE:-}" ]]; then
  cp -p "$SVCFIT_RUNTIME_LOCKFILE" "$provenance_dir/conda-explicit-linux-64.txt"
fi

{
  printf 'mode\tscope\tinput_dir\ttruth_dir\toutput_dir\texpected_cases\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$mode" "$scope" "$input_dir" "$truth_dir" "$output_dir" "$expected_cases"
} > "$provenance_dir/run-environment.tsv"
printf 'scenario\tpurity\tboot\texperiment\n' > "$provenance_dir/expected-cases.tsv"
for boot in "${boots[@]}"; do
  for purity in "${purities[@]}"; do
    for experiment in "${experiments[@]}"; do
      printf 'S1\t%s\t%s\t%s\n' "$purity" "$boot" "$experiment" >> "$provenance_dir/expected-cases.tsv"
    done
  done
done

slurm_site_args=()
[[ -n "${SLURM_ACCOUNT:-}" ]] && slurm_site_args+=(--account="$SLURM_ACCOUNT")
[[ -n "${SLURM_PARTITION:-}" ]] && slurm_site_args+=(--partition="$SLURM_PARTITION")
[[ -n "${SLURM_QOS:-}" ]] && slurm_site_args+=(--qos="$SLURM_QOS")

export HOME="$run_root/runtime-home"
export XDG_CACHE_HOME="$run_root/cache"
export PIP_CACHE_DIR="$run_root/cache/pip"
export R_USER="$run_root/r-user"
export R_LIBS_USER="$svcfit_r_lib"
export SCRIPT_DIR="$script_dir" INPUT_DIR="$input_dir" TRUTH_DIR="$truth_dir"
export OUTPUT_DIR="$output_dir" EVAL_DIR="$eval_dir" RUN_ROOT="$run_root"
export SVCFIT_REPO="$svcfit_repo" EXPECTED_SVCFIT_COMMIT="$svcfit_commit"
export ENV_SETUP_SCRIPT="$runtime_loader" RSCRIPT_BIN="${SVCFIT_R:?Set SVCFIT_R in the site configuration}"
export RETICULATE_PYTHON="${RETICULATE_PYTHON:?Set RETICULATE_PYTHON in the site configuration}"
export PYTHON_ENV= SCENARIO=S1 STAGES=svcfit+cluster+tree COV=50

job_ids=()
for boot in "${boots[@]}"; do
  for purity in "${purities[@]}"; do
    export BOOT="$boot" ppur="$purity"
    job=$("$sbatch_bin" --parsable "${slurm_site_args[@]}" \
      --chdir="$script_dir" --array="$array_spec" \
      --mem="${WORKER_MEM:-4G}" --time="${WORKER_TIME:-01:00:00}" \
      --output="$log_dir/longi_svcfit_%A_%a.out" \
      --error="$log_dir/longi_svcfit_%A_%a.err" --export=ALL \
      "$script_dir/longi_svcfit.sh")
    job=${job%%;*}
    job_ids+=("$job")
    printf 'worker\tboot=%s\tpurity=%s\tarray=%s\tjob=%s\n' "$boot" "$purity" "$array_spec" "$job" >> "$provenance_dir/submitted-jobs.tsv"
  done
done

dependency=$(IFS=:; printf 'afterok:%s' "${job_ids[*]}")
export EXPECTED_CASES="$expected_cases"
export EXPECTED_NO_TREE="$expected_no_tree"
export EXPECTED_CORRECT_TOPOLOGIES="$expected_correct"
export EVAL_PURITIES="$(IFS=,; printf '%s' "${purities[*]}")"
export EVAL_EXPERIMENTS="$(IFS=,; printf '%s' "${experiments[*]}")"
export N_BOOTS=5
if [[ "$scope" == smoke ]]; then export EVAL_BOOT=0; else export EVAL_BOOT=; fi

eval_job=$("$sbatch_bin" --parsable "${slurm_site_args[@]}" \
  --dependency="$dependency" --chdir="$script_dir" \
  --job-name="tree-eval-${mode}-${scope}" \
  --output="$log_dir/evaluate_%j.out" --error="$log_dir/evaluate_%j.err" \
  --mem="${EVAL_MEM:-6G}" --time="${EVAL_TIME:-00:30:00}" --export=ALL \
  "$script_dir/evaluate_downstream.sbatch")
eval_job=${eval_job%%;*}
printf 'evaluation\tdependency=%s\tjob=%s\n' "$dependency" "$eval_job" >> "$provenance_dir/submitted-jobs.tsv"

printf 'Evaluation job:  %s\n' "$eval_job"
printf 'Job record:      %s\n' "$provenance_dir/submitted-jobs.tsv"
