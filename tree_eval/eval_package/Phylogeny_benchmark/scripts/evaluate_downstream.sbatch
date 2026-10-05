#!/usr/bin/bash
#SBATCH --nodes=1

set -Eeuo pipefail

: "${SCRIPT_DIR:?}"
: "${INPUT_DIR:?}"
: "${TRUTH_DIR:?}"
: "${OUTPUT_DIR:?}"
: "${EVAL_DIR:?}"
: "${RUN_ROOT:?}"
: "${ENV_SETUP_SCRIPT:?}"
: "${RSCRIPT_BIN:?}"
: "${EXPECTED_CASES:?}"
: "${EVAL_PURITIES:?}"
: "${EVAL_EXPERIMENTS:?}"

# shellcheck disable=SC1090
source "$ENV_SETUP_SCRIPT"

args=(
  --work_dir "$OUTPUT_DIR"
  --input_dir "$INPUT_DIR"
  --truth_dir "$TRUTH_DIR"
  --out_dir "$EVAL_DIR"
  --n_boots "${N_BOOTS:-5}"
  --expected_cases "$EXPECTED_CASES"
  --purities "$EVAL_PURITIES"
  --experiments "$EVAL_EXPERIMENTS"
)
[[ -n "${EVAL_BOOT:-}" ]] && args+=(--boot "$EVAL_BOOT")
[[ -n "${EXPECTED_NO_TREE:-}" ]] && args+=(--expected_no_tree "$EXPECTED_NO_TREE")
[[ -n "${EXPECTED_CORRECT_TOPOLOGIES:-}" ]] && args+=(--expected_correct_topologies "$EXPECTED_CORRECT_TOPOLOGIES")

"$RSCRIPT_BIN" "$SCRIPT_DIR/evaluate_downstream.R" "${args[@]}"

cd "$RUN_ROOT"
find output evaluation -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > provenance/run-files.sha256
{
  printf 'completed_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  printf 'expected_cases\t%s\n' "$EXPECTED_CASES"
  printf 'evaluation_summary_sha256\t'
  sha256sum evaluation/evaluation_summary.csv | awk '{print $1}'
} > provenance/observed-results.tsv
