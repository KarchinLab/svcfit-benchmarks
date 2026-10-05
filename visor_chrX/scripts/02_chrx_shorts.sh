#!/bin/bash
#SBATCH --job-name=chrx_shorts
#SBATCH --output=log/shorts_%A_%a.out
#SBATCH --error=log/shorts_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=8:00:00
#
# NO --array DIRECTIVE HERE, DELIBERATELY. It used to be `--array=0-44%45`, which
# applies no matter what argument is passed — so `sbatch 02_chrx_shorts.sh normal`
# launched 45 concurrent copies of the *matched normal* job, all writing to the
# same output directory. 43 died on VISOR's "output folder is not empty" guard
# (SHORtS.py:880) and the two survivors clobbered each other (job 29344043, §28).
# The array size now lives on the command line, so the two modes cannot be
# confused. Pass --array ONLY for the tumour run.
#
# SIZING (see §26). The 45 tasks are fully independent, so concurrency is the
# real lever: %45 runs them in one wave instead of 5 sequential batches.
# 45 x 12 = 540 cpus against a `normal` QoS ceiling of cpu=3600, on 48-core nodes
# (4 tasks/node). Cores-per-task help LESS than they look: SHORtS.py:1054 walks
# regions serially ("do not use multi-processing on this as minimap2 may require
# too much memory") and wgsim read generation is single-threaded. --threads only
# reaches minimap2 -t, samtools sort -@ threads/2 and samtools merge -@. So the
# read-simulation half of each task does not speed up at all.
#
# Simulate short reads for the hemizygous chrX simulation.
#
#   sbatch scripts/02_chrx_shorts.sh normal                    <- matched normal, ONE job, run FIRST
#   sbatch --array=0-44%45 scripts/02_chrx_shorts.sh           <- 45 tumour BAMs (3 exps x 15 conds)
#
# The --array belongs on the command line, NOT in a directive: a directive would
# also apply to the `normal` invocation and fan it out 45 ways (§28).
#
# The matched normal is a separate single job because it is condition-independent:
# one germline BAM serves every tumour. Run it first — it is also the cleanest
# place to verify the chrX/chr22 = 0.5 hemizygosity signature, since it carries
# no SVs or CNVs to perturb the ratio.
#
# The most dangerous error available here: the matched
# normal must be single-copy on chrX too. cn_bar is a T/N depth ratio, so a
# diploid normal against a haploid tumour halves every cn_bar and every chrX
# fraction comes out 2x wrong with nothing appearing broken. Both the matched
# normal AND the normal-contamination clone inside the tumour BAM (§12a) are
# built from clone_genomes/normal/, which is haploid-chrX by construction (§22).

set -euo pipefail

# LOCATING chrx_common.sh UNDER SBATCH.
# sbatch copies this script to a node-local spool dir and runs it from there, so
# ${BASH_SOURCE[0]} is /cm/local/apps/slurm/var/spool/job<N>/slurm_script and a
# BASH_SOURCE-relative source() fails with "No such file or directory". That is
# exactly what killed all 45 tasks of job 29342365 (§27). The submitted
# 00_visor_shorts.sh:20 carries a SLURM_SUBMIT_DIR fallback for the same reason.
# Try, in order: explicit override, submit dir, BASH_SOURCE dir, install path.
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

# A batch job does not reliably inherit an activated Conda environment, so
# activate the configured absolute project prefix here.
if ! command -v VISOR >/dev/null 2>&1; then
    # CONDA_SH and ENV_VISOR come from config.local.sh via chrx_common.sh.
    [[ -n "${CONDA_SH:-}" ]] || { echo "ERROR: CONDA_SH unset; is config.local.sh loaded?" >&2; exit 1; }
    [[ -r "$CONDA_SH" ]] || { echo "ERROR: no VISOR on PATH and cannot read $CONDA_SH" >&2; exit 1; }
    set +u; source "$CONDA_SH"; conda activate "${VISOR_ENV:-${ENV_VISOR}}"; set -u
fi

command -v VISOR    >/dev/null || { echo "ERROR: VISOR not on PATH after activating ${VISOR_ENV:-visor}" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }
[[ -r "$SEEDED_PY" ]] || { echo "ERROR: cannot read $SEEDED_PY" >&2; exit 1; }

MODE="${1:-tumour}"

# VISOR emits non-fatal warnings; we inspect the log explicitly rather than
# letting a pipe status decide (§21).
set +o pipefail

# ---------------------------------------------------------------------------
run_shorts() {
    local outdir="$1" seed="$2" bed="$3" log="$4"; shift 4
    local sources=("$@")
    local fracs=("${FRACS[@]}")

    mkdir -p "$outdir"
    rm -rf "$outdir/short.out"

    echo "  sources : ${sources[*]}"
    echo "  fracs   : ${fracs[*]}"
    echo "  coverage: $cov   seed: $seed"

    set +e
    "$ENV_VISOR/bin/python" "$SEEDED_PY" "$seed" \
        -g "$REF" \
        -s "${sources[@]}" \
        -b "$bed" \
        -o "$outdir/short.out" \
        --coverage "$cov" \
        --clonefraction "${fracs[@]}" \
        --mutation 0 \
        --error 0 \
        --indels 0 \
        --extindels 0 \
        --threads "${SLURM_CPUS_PER_TASK:-2}" \
        > "$log" 2>&1
    local rc=$?
    set -e

    # VISOR EXIT CODE IS NOT TRUSTWORTHY ON /vast (§30).
    # After merging the final BAM, SHORtS os.remove()s the per-haplotype temp
    # files then os.rmdir()s their directories. On NFS the directory can still
    # report non-empty for a moment, so rmdir raises
    #   OSError: [Errno 39] Directory not empty: .../short.out/clone1/h1
    # and VISOR exits 1 — AFTER sim.srt.bam is complete and indexed. Job 29344957
    # died this way with a perfectly good 754 MB BAM.
    # So: tolerate that specific failure only when the BAM itself verifies.
    # Same lesson as §11 — check for the artefact, not the reported status.
    # User confirms this cleanup failure is not intermittent — VISOR fails it
    # every run on this filesystem. So match the condition broadly (errno number
    # OR the message text, since the numeric errno is platform-specific) while
    # still requiring the BAM to verify. Anything else propagates.
    if [[ $rc -ne 0 ]]; then
        if grep -qE 'Errno 39|Directory not empty' "$log" \
           && samtools quickcheck "$outdir/short.out/sim.srt.bam" 2>/dev/null; then
            echo "  NOTE: VISOR exited $rc on its post-merge cleanup (rmdir on a"
            echo "        non-empty haplotype dir). Expected on this filesystem;"
            echo "        sim.srt.bam passes samtools quickcheck — continuing."
        else
            echo "ERROR: VISOR SHORtS failed (exit $rc)" >&2
            tail -20 "$log" >&2
            return "$rc"
        fi
    fi
}

# ===========================================================================
# MATCHED NORMAL
# ===========================================================================
if [[ "$MODE" == "normal" ]]; then
    # Belt and braces after §28: even if an --array sneaks back in, refuse to run
    # the matched normal as a fan-out. Every task would target one output dir,
    # tripping VISOR's non-empty-output guard or racing to corrupt the BAM.
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "ERROR: 'normal' is a single job — do not submit it with --array." >&2
        echo "       Run: sbatch ${BASH_SOURCE[0]} normal" >&2
        exit 1
    fi
    # COV-NAMED, like every tumour condition. It was "${short_base}/normal", which is
    # coverage-blind, while run_shorts() does `rm -rf "$outdir/short.out"` -- so simulating a
    # second coverage silently DELETED the first one's matched normal and left its committed
    # cn_bar and scores resting on a normal that no longer existed at that depth.
    outdir="${short_base}/normal_c${cov}"
    log="${outdir}/shorts.log"
    mkdir -p "$outdir"

    echo "=== matched normal BAM ==="
    build_short_bed "$outdir/short.bed" "$clone_genomes/normal"/*.fa.fai
    cat "$outdir/short.bed"

    FRACS=(100)
    run_shorts "$outdir" 999 "$outdir/short.bed" "$log" "$clone_genomes/normal"

    bam="$outdir/short.out/sim.srt.bam"
    [[ -s "$bam" ]] || { echo "ERROR: $bam not produced" >&2; tail -20 "$log" >&2; exit 1; }
    samtools index -@ "${SLURM_CPUS_PER_TASK:-2}" "$bam" 2>/dev/null || true

    echo "=== guards ==="
    # one clone dir, whose h2.fa deliberately lacks chrX -> exactly 1 warning
    check_warnings "$log" 1 || exit 1

    read -r ratio adepth xdepth < <(depth_ratio "$bam")
    printf '  %s depth %s | %s depth %s | ratio %s (expect ~0.5)\n' \
           "$AUTOSOME" "$adepth" "$HEMI" "$xdepth" "$ratio"
    awk -v r="$ratio" 'BEGIN { exit !(r > 0.47 && r < 0.53) }' || {
        echo "  GUARD FAILED: normal ${HEMI}/${AUTOSOME} depth ratio $ratio outside [0.47,0.53]." >&2
        echo "  A diploid normal here halves every cn_bar." >&2
        exit 1; }

    echo "RESULT: matched normal built and verified — $bam"
    exit 0
fi

# ===========================================================================
# TUMOUR ARRAY
# ===========================================================================
: "${SLURM_ARRAY_TASK_ID:?run under sbatch --array, or pass 'normal'}"
decode_task "$SLURM_ARRAY_TASK_ID"

echo "========================================================"
echo "  task      : ${SLURM_ARRAY_TASK_ID}   job: ${SLURM_JOB_ID:-NA}"
echo "  host      : $(hostname -s)"
echo "  experiment: ${exp_name}   condition: ${cond_name}"
echo "  purity    : ${ppur}%   mixture: ${sub_mix}%"
echo "  clone frac: c2=${frac[0]}% c3=${frac[1]}% normal=${frac[2]}%"
echo "  seed      : ${visor_seed}"
echo "========================================================"

# IDEMPOTENT, like 06 and 12. This stage had no skip guard, so resubmitting an
# array to recover a handful of tasks re-simulated all 45 -- ~28 task-hours to
# replace four BAMs. It matters more than convenience at 30 replicates: this
# cluster drops whole nodes' worth of tasks at start-up (exit 0:53, seconds in,
# no logs written), which was seen on three separate occasions during this run,
# and each recovery would otherwise cost a full replicate's simulation.
#
# The guard tests the INDEX, not just the BAM: an interrupted task leaves a
# partial sim.srt.bam behind, and `rm -rf short.out` at the top of run_shorts
# means a BAM without its .bai is exactly the state a killed task leaves. The
# index is written last, so its presence is what says the task finished.
bam="${work_dir}/short.out/sim.srt.bam"
if [[ -s "$bam" && -s "${bam}.bai" ]]; then
    echo "SKIP ${exp_name}/${cond_name}: $bam already present and indexed"
    exit 0
fi

c2dir="${clone_genomes}/${exp_name}/c2"
c3dir="${clone_genomes}/${exp_name}/c3"
ndir="${clone_genomes}/normal"
for d in "$c2dir" "$c3dir" "$ndir"; do
    [[ -s "$d/h1.fa.fai" && -s "$d/h2.fa.fai" ]] || { echo "ERROR: missing genomes in $d" >&2; exit 1; }
done

mkdir -p "$work_dir"
log="${work_dir}/shorts.log"
build_short_bed "$work_dir/short.bed" "$c2dir"/*.fa.fai "$c3dir"/*.fa.fai "$ndir"/*.fa.fai

FRACS=("${frac[@]}")
run_shorts "$work_dir" "$visor_seed" "$work_dir/short.bed" "$log" "$c2dir" "$c3dir" "$ndir"

[[ -s "$bam" ]] || { echo "ERROR: $bam not produced" >&2; tail -20 "$log" >&2; exit 1; }

# Indexing was `|| true`. It cannot be, now that the skip guard above treats the
# .bai as the marker of a finished task: a swallowed failure would leave a task
# that exits 0, is re-run by every later sweep, and hands Manta and `samtools
# depth` a BAM they cannot open. Both need the index regardless.
samtools index -@ "${SLURM_CPUS_PER_TASK:-2}" "$bam" || {
    echo "ERROR: samtools index failed on $bam" >&2; exit 1; }
[[ -s "${bam}.bai" ]] || { echo "ERROR: no ${bam}.bai after indexing" >&2; exit 1; }

echo "=== guards ==="
# three clone dirs, each with an h2.fa that deliberately lacks chrX
check_warnings "$log" 3 || exit 1

# Reported, not asserted: the tumour ratio is legitimately perturbed by the SVs
# and (in e2/e4) by large amplifications, so there is no single expected value.
# cn_bar is estimated downstream from depth (05_chrx_depth_segmentation_sim.R).
read -r ratio adepth xdepth < <(depth_ratio "$bam")
printf '  %s depth %s | %s depth %s | ratio %s (informational)\n' \
       "$AUTOSOME" "$adepth" "$HEMI" "$xdepth" "$ratio"

echo "RESULT: ${exp_name}/${cond_name} — $bam ($(du -sh "$bam" | cut -f1))"
