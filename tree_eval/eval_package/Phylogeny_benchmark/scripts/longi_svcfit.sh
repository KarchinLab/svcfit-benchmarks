#!/usr/bin/env bash
#SBATCH --job-name=longi_svcfit
#SBATCH --output=log/longi_svcfit_%A_%a.out
#SBATCH --error=log/longi_svcfit_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --time=1:00:00
#SBATCH --array=0-4

set -Eeuo pipefail

###############################################################################
# Run SVCFit, clustering, and tree reconstruction for one experiment in an
# S1 purity/bootstrap group. The driver exports SCRIPT_DIR because SLURM runs a
# spooled copy of this file, making BASH_SOURCE point at the spool directory.
###############################################################################

cov=${COV:-50}
SCENARIO=${SCENARIO:-S1}
ppur=${ppur:-60}
BOOT=${BOOT:-0}
STAGES=${STAGES:-svcfit+cluster+tree}
STAGES=${STAGES//+/,}

: "${SCRIPT_DIR:?SCRIPT_DIR must be exported by the submission driver}"
: "${INPUT_DIR:?INPUT_DIR must point to the existing simulation/calling inputs}"
: "${OUTPUT_DIR:?OUTPUT_DIR must point to a new run output directory}"
: "${SVCFIT_REPO:?SVCFIT_REPO must point to the pinned SVCFit checkout}"
: "${EXPECTED_SVCFIT_COMMIT:?EXPECTED_SVCFIT_COMMIT is required}"
RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}
ENV_SETUP_SCRIPT=${ENV_SETUP_SCRIPT:-}
PYTHON_ENV=${PYTHON_ENV-py3}

[[ -f "$SCRIPT_DIR/longi_svcfit.R" ]] || {
  echo "ERROR: longi_svcfit.R not found at $SCRIPT_DIR/longi_svcfit.R" >&2
  exit 1
}
[[ -d "$INPUT_DIR" ]] || {
  echo "ERROR: INPUT_DIR is not a directory: $INPUT_DIR" >&2
  exit 1
}
git -C "$SVCFIT_REPO" rev-parse --is-inside-work-tree >/dev/null 2>&1 || {
  echo "ERROR: SVCFIT_REPO is not a Git checkout or worktree: $SVCFIT_REPO" >&2
  exit 1
}

actual_svcfit_commit=$(git -C "$SVCFIT_REPO" rev-parse HEAD)
[[ "$actual_svcfit_commit" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
  echo "ERROR: SVCFit commit is $actual_svcfit_commit; expected $EXPECTED_SVCFIT_COMMIT" >&2
  exit 1
}

exp_lst=(exp1 exp2 exp3 exp4 exp5)
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
[[ "$task_id" =~ ^[0-4]$ ]] || {
  echo "ERROR: SLURM_ARRAY_TASK_ID must be 0-4: $task_id" >&2
  exit 1
}
samp_name=${exp_lst[$task_id]}

case_rel="$SCENARIO/c${cov}p${ppur}"
if [[ "$BOOT" -eq 0 ]]; then
  case_rel="$case_rel/$samp_name"
else
  case_rel="$case_rel/b${BOOT}/$samp_name"
fi
case_output="$OUTPUT_DIR/$case_rel"
mkdir -p "$case_output"
run_stamp="$case_output/.job-start-${SLURM_JOB_ID:-manual}-${task_id}"
: > "$run_stamp"

if [[ -n "$ENV_SETUP_SCRIPT" ]]; then
  [[ -r "$ENV_SETUP_SCRIPT" ]] || {
    echo "ERROR: ENV_SETUP_SCRIPT is not readable: $ENV_SETUP_SCRIPT" >&2
    exit 1
  }
  # shellcheck source=/dev/null
  source "$ENV_SETUP_SCRIPT"
fi
command -v "$RSCRIPT_BIN" >/dev/null 2>&1 || {
  echo "ERROR: Rscript command is unavailable after environment setup: $RSCRIPT_BIN" >&2
  exit 1
}

echo "=== Pipeline start: $SCENARIO / $samp_name / purity=$ppur / boot=$BOOT ==="
echo "Stages:        $STAGES"
echo "Input root:    $INPUT_DIR"
echo "Output root:   $OUTPUT_DIR"
echo "SVCFit commit: $actual_svcfit_commit"

"$RSCRIPT_BIN" "$SCRIPT_DIR/longi_svcfit.R" \
  --scenario "$SCENARIO" \
  --purity "$ppur" \
  --exp "$samp_name" \
  --input_dir "$INPUT_DIR" \
  --output_dir "$OUTPUT_DIR" \
  --svcfit_repo "$SVCFIT_REPO" \
  --python_env "$PYTHON_ENV" \
  --cov "$cov" \
  --stages "$STAGES" \
  --boot "$BOOT"

require_fresh() {
  local path=$1
  [[ -s "$path" ]] || {
    echo "ERROR: required output is absent or empty: $path" >&2
    exit 1
  }
  [[ "$path" -nt "$run_stamp" ]] || {
    echo "ERROR: output is not newer than the job start marker: $path" >&2
    exit 1
  }
}

if [[ ",$STAGES," == *,svcfit,* ]]; then
  require_fresh "$case_output/svcfit_output/${SCENARIO}_p${ppur}_${samp_name}_t1.bed"
  require_fresh "$case_output/svcfit_output/${SCENARIO}_p${ppur}_${samp_name}_t2.bed"
fi

if [[ ",$STAGES," == *,cluster,* || ",$STAGES," == *,tree,* ]]; then
  require_fresh "$case_output/clustering/sv2cluster.csv"
  require_fresh "$case_output/clustering/cluster_result.csv"
  require_fresh "$case_output/clustering/cluster_centroids.csv"
fi

if [[ ",$STAGES," == *,tree,* ]]; then
  tree_file="$case_output/tree/tree_result.rds"
  no_tree_file="$case_output/tree/NO_TREE.txt"
  if [[ -s "$tree_file" && -s "$no_tree_file" ]]; then
    echo "ERROR: both tree_result.rds and NO_TREE.txt exist: $case_output/tree" >&2
    exit 1
  fi
  if [[ -s "$tree_file" ]]; then
    require_fresh "$tree_file"
    require_fresh "$case_output/tree/tree_edges.csv"
    require_fresh "$case_output/tree/topology.txt"
  elif [[ -s "$no_tree_file" ]]; then
    require_fresh "$no_tree_file"
  else
    echo "ERROR: tree stage produced neither a tree nor NO_TREE.txt: $case_output/tree" >&2
    exit 1
  fi
fi

mv "$run_stamp" "$case_output/.job-complete-${SLURM_JOB_ID:-manual}-${task_id}"
echo "=== Done: $SCENARIO / $samp_name / purity=$ppur / boot=$BOOT ==="
