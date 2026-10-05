#!/bin/bash
#SBATCH --job-name=chrx_call
#SBATCH --output=log/call_%A_%a.out
#SBATCH --error=log/call_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#
# Manta + SVtyper across the condition grid, one task per tumour BAM.
# Wraps 11_chrx_manta_svtyper.sh; see that script for why each step matches
# visor_replicates/01_manta_svtyper.sh.
#
# NO --array DIRECTIVE HERE, DELIBERATELY -- same reason as 06:
#
#   sbatch --array=0-44 scripts/12_chrx_calls_array.sbatch    # all 45 (idempotent)
#
# --time is 2:00:00 = TWO HOURS, long-form on purpose (gotcha 5). Observed pilot
# runtime is ~40 s/task, so this is very generous; Manta's own retry logic can
# stretch a task if a workflow step fails.
#
# IDEMPOTENT: a task whose genotyped VCF already exists and is non-empty exits 0
# without recalling, so a partially-failed array can be resubmitted wholesale.
# Manta additionally refuses to reuse a populated runDir, so the manta directory
# is cleared before a real (non-skipped) run.

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
             "${SLURM_SUBMIT_DIR:-}/scripts" \
             "${SLURM_SUBMIT_DIR:-}" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" \
             "$(_visor_walk_up)"; do
        [[ -n "$d" && -r "$d/11_chrx_manta_svtyper.sh" ]] && { printf '%s' "$d"; return 0; }
    done
    return 1
}
SCRIPTS="$(_resolve_scripts_dir)" || { echo "ERROR: cannot locate scripts dir" >&2; exit 1; }

# For cov and cond_tsv. The conditions file is coverage-named, so the coverage
# this array runs is fixed by $COV at submit time and cannot disagree with the
# normal BAM that 11 pairs against:
#   sbatch --export=ALL,COV=50 --array=0-44%45 scripts/12_chrx_calls_array.sbatch
source "${SCRIPTS}/chrx_common.sh"

conds="$cond_tsv"
[[ -s "$conds" ]] || { echo "ERROR: missing $conds (run COV=${cov} 04_chrx_ground_truth.sh)" >&2; exit 1; }

tid="${SLURM_ARRAY_TASK_ID:?this script must run as an array task}"
read -r _tid exp cond _rest < <(awk -F'\t' -v t="$tid" 'NR>1 && $1==t' "$conds")
[[ -n "${exp:-}" && -n "${cond:-}" ]] || { echo "ERROR: no condition for task $tid" >&2; exit 1; }

name="${exp}_${cond}"
# $CALLS_DIR, not ${BASE_DIR}/calls: chrx_common.sh rep-tags it, and 11 writes
# there. Hardcoding the untagged path meant that under REP>0 the skip test below
# consulted REPLICATE 0's VCFs, which exist -- so every condition would report
# SKIP and the replicate would end with no calls of its own, while `rm -rf
# "${out}/manta"` deleted replicate 0's Manta run directories on the way past.
out="${CALLS_DIR}/${name}"
vcf="${out}/svtyp/svt_${name}.vcf"

if [[ -s "$vcf" ]] && grep -qv '^#' "$vcf" 2>/dev/null; then
    echo "SKIP ${name}: $vcf already present and non-empty"
    exit 0
fi

rm -rf "${out}/manta"          # Manta will not reuse a populated runDir
echo "task ${tid} -> ${name}   ($(date))"

# NOT `exec`, so this can check its own post-condition: a task that produced no
# genotyped VCF must fail.
"${SCRIPTS}/11_chrx_manta_svtyper.sh" "$exp" "$cond" "${SLURM_CPUS_PER_TASK:-4}"
rc=$?

[[ $rc -eq 0 ]] || { echo "ERROR ${name}: 11_chrx_manta_svtyper.sh exited $rc" >&2; exit "$rc"; }
[[ -s "$vcf" ]] || { echo "ERROR ${name}: exited 0 but produced no $vcf" >&2; exit 1; }
echo "OK ${name}: $(grep -vc '^#' "$vcf" || true) genotyped records"
