#!/bin/bash
# Reconcile replicates against what is actually on disk, and resubmit the gaps.
#
#   COV=50 ./18_chrx_replicate_sweep.sh [FIRST] [LAST] [--dry-run]
#
# WHY THIS EXISTS. The chains in 15 wire their stages with `afterok`, which is
# correct in principle -- there is nothing to call or segment if the reads did
# not simulate -- but it is all-or-nothing across a 45-task array. This cluster
# drops whole nodes' worth of tasks at start-up (exit 0:53, a few seconds in, no
# log files written at all, seen on three separate occasions during this run),
# and four lost tasks out of 45 are enough to leave an entire replicate's
# downstream permanently PENDING on an afterok that can never be satisfied.
#
# So completeness is decided by looking at the ARTEFACTS rather than by trusting
# job states. Every stage is idempotent, so this is safe to run repeatedly. It
# resubmits only the missing array indices, which is the difference between
# re-simulating four BAMs and re-simulating all forty-five.
#
# Run it while jobs are still in flight and it will re-submit work that is
# already queued. That is wasteful but not wrong -- the stages skip completed
# outputs -- so prefer to sweep when the queue for those replicates is empty.

set -euo pipefail

DRY=""
args=()
for a in "$@"; do [[ "$a" == "--dry-run" ]] && DRY=1 || args+=("$a"); done
FIRST="${args[0]:-1}"; LAST="${args[1]:-30}"

COV="${COV:?set COV, e.g. COV=50 $0}"
_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load config before selecting the data root; chrx_common.sh defaults BASE_DIR
# to CHRX_DIR rather than the source checkout.
source "${_DIR}/chrx_common.sh"
cd "$BASE_DIR"

if ! command -v sbatch >/dev/null 2>&1 && [[ -n "${SLURM_BIN:-}" ]]; then
    export PATH="$SLURM_BIN:$PATH"
fi
command -v sbatch >/dev/null 2>&1 || {
    echo "ERROR: sbatch is not on PATH; set SLURM_BIN in config.local.sh" >&2
    exit 1
}
conds="${BASE_DIR}/truth/chrx_conditions_c${COV}.tsv"
[[ -s "$conds" ]] || { echo "ERROR: missing $conds" >&2; exit 1; }

# For CHRX_EXCLUDE_NODES -- the nodes that kill tasks at start-up. Defined in
# chrx_common.sh rather than here so 15 and 18 cannot drift apart on it. The
# sweeper is the script that most needs it: it exists to repair exactly the
# damage those nodes do, so submitting the repairs back onto them is circular.
excl=""; [[ -n "${CHRX_EXCLUDE_NODES:-}" ]] && excl="--exclude=${CHRX_EXCLUDE_NODES}"

# REFUSE TO SWEEP OVER A LIVE QUEUE. Completeness here is judged from artefacts,
# and a replicate that is merely still running is indistinguishable from one that
# failed -- both show missing BAMs. Sweeping mid-flight therefore resubmits work
# that is already queued: harmless to the results, since every stage is
# idempotent, but it doubles a 540-cpu-per-replicate simulation and makes the
# queue unreadable.
#
# --force is for the deliberate case where the remaining jobs are known to be
# stranded on an afterok that cannot be satisfied.
live=$(squeue -u "$(whoami)" -h -o "%j" 2>/dev/null | grep -c '^chrx_' || true)
if (( live > 0 )) && [[ -z "$DRY" && -z "${FORCE:-}" ]]; then
    echo "REFUSING: ${live} chrx_* job(s) still in the queue." >&2
    echo "  A replicate that is still running looks exactly like one that failed here," >&2
    echo "  so sweeping now would resubmit work already in flight." >&2
    echo "  Wait for the queue to drain, or re-run with FORCE=1 if those jobs are" >&2
    echo "  known to be stranded (e.g. PENDING on DependencyNeverSatisfied)." >&2
    exit 1
fi

printf '%-4s %-8s %-8s %-8s %-6s  %s\n' rep BAMs VCFs segs score action
printf '%-4s %-8s %-8s %-8s %-6s  %s\n' ---- -------- -------- -------- ------ ------

resubmitted=0
for (( rep = FIRST; rep <= LAST; rep++ )); do
    rt=""; (( rep > 0 )) && rt="_rep${rep}"
    short_d="${BASE_DIR}/short${rt}"
    calls_d="${BASE_DIR}/calls${rt}"
    seg_d="${BASE_DIR}/seg10k${rt}"
    score_f="${BASE_DIR}/scoring${rt}/chrx_svcfit_scores_c${COV}.tsv"

    miss_bam=(); miss_vcf=(); miss_seg=()
    while IFS=$'\t' read -r tid exp cond _rest; do
        b="${short_d}/${exp}/${cond}/short.out/sim.srt.bam"
        [[ -s "$b" && -s "${b}.bai" ]]                                        || miss_bam+=("$tid")
        [[ -s "${calls_d}/${exp}_${cond}/svtyp/svt_${exp}_${cond}.vcf" ]]     || miss_vcf+=("$tid")
        [[ -s "${seg_d}/chrx_seg_${exp}_${cond}.csv" ]]                       || miss_seg+=("$tid")
    done < <(tail -n +2 "$conds")

    nb=${#miss_bam[@]}; nv=${#miss_vcf[@]}; ns=${#miss_seg[@]}
    have_score=$([[ -s "$score_f" ]] && echo yes || echo no)

    if (( nb == 0 && nv == 0 && ns == 0 )) && [[ "$have_score" == yes ]]; then
        action="complete"
    elif (( nb == 0 && nv == 0 && ns == 0 )); then
        action="score only"
    else
        action="resubmit ${nb}/${nv}/${ns}"
    fi

    printf '%-4s %-8s %-8s %-8s %-6s  %s\n' \
        "$rep" "$((45-nb))/45" "$((45-nv))/45" "$((45-ns))/45" "$have_score" "$action"

    [[ "$action" == "complete" ]] && continue
    [[ -n "$DRY" ]] && continue
    resubmitted=$(( resubmitted + 1 ))

    EXP="ALL,COV=${COV},REP=${rep}"
    mkdir -p "$short_d"
    ln -sfn "${BASE_DIR}/short/normal_c${COV}" "${short_d}/normal_c${COV}"

    # afterany, not afterok, throughout: a stage that lost tasks must not strand
    # the stages behind it a second time. Each one asserts its own inputs and
    # fails per task, which a later sweep then picks up.
    dep=""; wait_on=()
    if (( nb > 0 )); then
        list=$(IFS=,; echo "${miss_bam[*]}")
        a=$(sbatch --parsable $excl --export="$EXP" --array="${list}%45" scripts/02_chrx_shorts.sh)
        echo "       shorts ${a}  <- ${nb} task(s): ${list}"
        dep="--dependency=afterany:${a}"
    fi
    if (( nv > 0 )); then
        list=$(IFS=,; echo "${miss_vcf[*]}")
        b=$(sbatch --parsable $dep $excl --export="$EXP" --array="${list}%45" \
                scripts/12_chrx_calls_array.sbatch)
        echo "       calls  ${b}  <- ${nv} task(s)"
        wait_on+=("$b")
    fi
    if (( ns > 0 )); then
        list=$(IFS=,; echo "${miss_seg[*]}")
        c=$(sbatch --parsable $dep $excl --export="${EXP},BIN_KB=10,SEG_DIR=${BASE_DIR}/seg10k" \
                --time=6:00:00 --array="${list}%45" scripts/06_chrx_segment_array.sbatch)
        echo "       seg    ${c}  <- ${ns} task(s)"
        wait_on+=("$c")
    fi

    # Scoring needs a COMPLETE replicate, so it is the one place afterok is
    # right: a partial input here would silently score fewer conditions than the
    # replicate is supposed to have, and 09 skips a condition whose segmentation
    # is missing rather than failing on it.
    if [[ "$have_score" == no ]]; then
        sdep=""
        (( ${#wait_on[@]} )) && sdep="--dependency=afterok:$(IFS=:; echo "${wait_on[*]}")"
        d=$(sbatch --parsable $sdep $excl --export="$EXP" scripts/17_chrx_score_one.sbatch)
        echo "       score  ${d}"
    fi
done

echo
if [[ -n "$DRY" ]]; then
    echo "dry run: nothing submitted."
else
    echo "${resubmitted} replicate(s) touched. Re-run after these land; it is idempotent."
fi
