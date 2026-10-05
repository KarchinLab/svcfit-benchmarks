#!/bin/bash
#SBATCH --job-name=chrx_seg
#SBATCH --output=log/seg_%A_%a.out
#SBATCH --error=log/seg_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#
# chrX depth segmentation across the condition grid, one task per tumour BAM.
#
# NO --array DIRECTIVE HERE, DELIBERATELY. A directive applies
# regardless of arguments, which is how the matched-normal job once got fanned out
# 45 ways onto one output directory. Array membership belongs on the command line:
#
#   sbatch --array=0-11,13-41,43-44 scripts/06_chrx_segment_array.sbatch   # the 43 not yet done
#   sbatch --array=0-44             scripts/06_chrx_segment_array.sbatch   # all 45 (idempotent)
#
# --time above is 1:00:00 = ONE HOUR, written long-form on purpose: `--time=120:00`
# would be 120 MINUTES, not hours.
# Observed runtime is ~2 min/task, so an hour is generous.
#
# Each task is independent and writes its own CSV, so there is no shared output to
# race on. Tasks re-read the same matched-normal BAM; that is redundant work the
# cohort script also does, left alone rather than optimised so this stays diffable
# against prostate_mixture/scripts/helper/chrx_depth_segmentation.R.
#
# IDEMPOTENT: a task whose CSV already exists and is non-empty exits 0 without
# recomputing, so a partially-failed array can be resubmitted wholesale.

set -euo pipefail

# sbatch copies this script to a node-local spool dir, so ${BASH_SOURCE[0]} is
# /cm/local/.../slurm_script and a BASH_SOURCE-relative source() fails (gotcha 2).
# Last-resort locator, replacing a hardcoded
# a checkout-specific scripts path that only ever existed on
# one of the three machines. Walks up from the submit directory to the repo root
# (identified by config.local.sh) and derives the scripts dir from there, so it
# works wherever the checkout lives. BASH_SOURCE is useless under sbatch, which
# copies the script to a node-local spool dir.
_visor_walk_up() {
    local d="${SLURM_SUBMIT_DIR:-$PWD}"
    while [[ "$d" != "/" && -n "$d" ]]; do
        [[ -r "$d/config.local.sh" ]] && { printf '%s' "$d/visor_chrX/scripts"; return 0; }
        d="$(dirname "$d")"
    done
    return 1
}

_resolve_scripts_dir() {
    local d
    for d in "${CHRX_SCRIPTS:-}" \
             "${SLURM_SUBMIT_DIR:-}" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" \
             "$(_visor_walk_up)"; do
        [[ -n "$d" && -r "$d/chrx_common.sh" ]] && { printf '%s' "$d"; return 0; }
    done
    return 1
}
_DIR="$(_resolve_scripts_dir)" || {
    echo "ERROR: cannot locate chrx_common.sh — set CHRX_SCRIPTS=/path/to/scripts" >&2; exit 1; }
source "${_DIR}/chrx_common.sh"

# A batch job does not reliably inherit an activated Conda environment, so use
# the configured project prefix. conda.sh trips `set -u`.
if ! command -v Rscript >/dev/null 2>&1 || ! command -v samtools >/dev/null 2>&1; then
    # CONDA_SH and ENV_VISOR come from config.local.sh via chrx_common.sh.
    [[ -n "${CONDA_SH:-}" ]] || { echo "ERROR: CONDA_SH unset; is config.local.sh loaded?" >&2; exit 1; }
    [[ -r "$CONDA_SH" ]] || { echo "ERROR: cannot read $CONDA_SH" >&2; exit 1; }
    set +u; source "$CONDA_SH"; conda activate "$ENV_VISOR"; set -u
fi
command -v Rscript  >/dev/null || { echo "ERROR: Rscript not on PATH" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }

: "${SLURM_ARRAY_TASK_ID:?run under sbatch --array}"
decode_task "$SLURM_ARRAY_TASK_ID"

# Bin width and output directory are overridable so a finer re-segmentation does
# not clobber the 100 kb results. Both default to the original run.
#   BIN_KB=10 SEG_DIR=$V/seg10k sbatch --time=4:00:00 --array=0-44 06_...sbatch
# Runtime scales inversely with bin width: 05_..._sim.R spawns two `samtools
# depth` calls PER BIN, so 10 kb is 10x the subprocesses and ~10x the wall clock
# (~2 min -> ~20 min per task). Raise --time on the command line to match.
BIN_KB="${BIN_KB:-100}"
seg_dir="${SEG_DIR:-${BASE_DIR}/seg}${rep_tag}"   # rep-safe; see chrx_common.sh
mkdir -p "$seg_dir"

tbam="${work_dir}/short.out/sim.srt.bam"
nbam="${short_base}/normal_c${cov}/short.out/sim.srt.bam"   # cov-named; see 02_chrx_shorts.sh
out="${seg_dir}/chrx_seg_${exp_name}_${cond_name}.csv"

# psi_sample = 2 EXACTLY, not an estimate. chr22 is copy-neutral diploid in every
# clone by construction, so the sample's mean autosomal
# copy number is 2 by definition. On the real cohort this comes from FACETS and
# psi_sample/2 runs 0.87-1.72, which is why 05_..._sim.R applies it at all.
PSI=2

echo "========================================================"
echo "  task      : ${SLURM_ARRAY_TASK_ID}   job: ${SLURM_JOB_ID:-NA}"
echo "  host      : $(hostname -s)"
echo "  experiment: ${exp_name}   condition: ${cond_name}"
echo "  clone frac: c2=${frac[0]}% c3=${frac[1]}% normal=${frac[2]}%"
echo "  psi_sample: ${PSI}   bin: ${BIN_KB} kb   out dir: ${seg_dir}"
echo "========================================================"

if [[ -s "$out" ]]; then
    echo "SKIP: $out already exists and is non-empty"
    exit 0
fi

[[ -s "$tbam" ]] || { echo "ERROR: missing tumour BAM $tbam" >&2; exit 1; }
[[ -s "$nbam" ]] || { echo "ERROR: missing normal BAM $nbam" >&2; exit 1; }

Rscript "${_DIR}/05_chrx_depth_segmentation_sim.R" "$tbam" "$nbam" "$PSI" "$out" "$BIN_KB"

[[ -s "$out" ]] || { echo "ERROR: $out not produced" >&2; exit 1; }
echo "RESULT: ${exp_name}/${cond_name} -> $out ($(wc -l < "$out") lines incl. header)"
