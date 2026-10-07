#!/bin/bash
#SBATCH --job-name=chrx_score
#SBATCH --output=log/score_%j.out
#SBATCH --error=log/score_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#
# Join + score for ONE replicate: 09 (segment -> SV join, maximum overlap) then
# 13 (score against the truth table). Submitted by 15_chrx_replicate_chain.sh
# once that replicate's calling and segmentation arrays are both in.
#
# Needs COV and REP in the environment; both come through --export.

set -euo pipefail

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
             "${SLURM_SUBMIT_DIR:-}/scripts" \
             "${SLURM_SUBMIT_DIR:-}" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" \
             "$(_visor_walk_up)"; do
        [[ -n "$d" && -r "$d/chrx_common.sh" ]] && { printf '%s' "$d"; return 0; }
    done
    return 1
}
SCRIPTS="$(_resolve_scripts_dir)" || { echo "ERROR: cannot locate scripts dir" >&2; exit 1; }
source "${SCRIPTS}/chrx_common.sh"          # COV / REP -> cov, rep_tag, paths

V="$BASE_DIR"
echo "=== replicate ${REP}, coverage ${cov} ==="

# 09 and 13 both append rep_tag to these themselves, so pass the untagged names.
export COV REP
SEG_DIR="${V}/seg10k" OUT_DIR="${V}/cn_bar10k_m" "${SCRIPTS}/09_chrx_join_all.sh"

# SVCFit comes from SVCFIT_R_LIB (an SVCFit 7f32d81 library); 13 asserts the hemizygous entry points.
: "${SVCFIT_R_LIB:?set SVCFIT_R_LIB to an SVCFit library (see tools/run_svcfit_self_evaluation.sh)}"
export R_LIBS="${SVCFIT_R_LIB}:${RLIB}"
export CHRX_DIR="$V"

rscript="${SVCFIT_R:-$RSCRIPT}"
[[ -x "$rscript" ]] || { echo "ERROR: Rscript not executable at $rscript" >&2; exit 1; }
"$rscript" "${SCRIPTS}/13_chrx_score_svcfit.R"

out="${V}/scoring${rep_tag}/chrx_svcfit_scores_c${cov}.tsv"
[[ -s "$out" ]] || { echo "ERROR: replicate ${REP} produced no $out" >&2; exit 1; }
echo "RESULT: replicate ${REP} scored -> $out ($(( $(wc -l < "$out") - 1 )) rows)"
