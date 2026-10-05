#!/bin/bash
# Emit the condition-dependent ground truth for the hemizygous chrX simulation.
#
# This is the last piece of the stated deliverable. The BAMs (45 tumour + 1
# matched normal) and the SV -> clone assignment already exist; what was missing
# is SV -> FRACTION, which is condition-dependent.
#
# THE CONDITION GRID IS NOT RE-ENCODED HERE. This script sources chrx_common.sh
# and calls decode_task() — the same function that drove 02_chrx_shorts.sh — so
# the truth cannot drift from what was actually simulated. Re-typing the purity
# and mixture vectors into a second script is precisely how a truth table ends up
# silently describing a run that never happened.
#
# SCOPE: cellular fractions only. Copy number at each SV locus is deliberately
# NOT emitted — cn_bar is estimated from depth by 05_chrx_depth_segmentation_sim.R at
# scoring time.
#
# Usage:  ./04_chrx_ground_truth.sh          (interactive; needs no conda env)

set -euo pipefail

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

truth_dir="${BASE_DIR}/truth"
# cond_tsv and truth_tsv are coverage-named and come from chrx_common.sh, so this
# cannot overwrite another coverage's truth table. Run as:  COV=50 ./04_...sh
out_tsv="$truth_tsv"

command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 1. Condition table, straight out of decode_task().
#    45 tasks = 3 experiments x 15 conditions, matching the submitted array.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] emitting condition grid from decode_task()"
{
    printf 'task_id\texperiment\tcondition\tpurity_pct\tmixture_pct\tc2_pct\tc3_pct\tnormal_pct\tseed\tbam\n'
    for tid in $(seq 0 44); do
        decode_task "$tid"
        printf '%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n' \
            "$tid" "$exp_name" "$cond_name" "$ppur" "$sub_mix" \
            "${frac[0]}" "${frac[1]}" "${frac[2]}" "$visor_seed" \
            "${work_dir}/short.out/sim.srt.bam"
    done
} > "$cond_tsv"

echo "  wrote $cond_tsv ($(( $(wc -l < "$cond_tsv") - 1 )) conditions)"

# ---------------------------------------------------------------------------
# 2. Join against the SV -> clone assignment.
# ---------------------------------------------------------------------------
python3 "${_DIR}/04_chrx_ground_truth.py" \
    --conditions   "$cond_tsv" \
    --assignment   "${truth_dir}/sv_beds/sv_subclone_assignment.tsv" \
    --sv-bed-dir   "${truth_dir}/sv_beds" \
    --out          "$out_tsv"
