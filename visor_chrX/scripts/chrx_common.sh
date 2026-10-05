#!/bin/bash
# Shared configuration for the hemizygous chrX simulation (session log §20-§23).
#
# 3 experiments x 15 conditions = 45 array tasks (0-44).
#   exp_idx  = SLURM_ARRAY_TASK_ID / 15    (0-2 -> e1, e2, e4)
#   cond_idx = SLURM_ARRAY_TASK_ID % 15    (0-14)
#
# Condition index -> (purity, mixture), same grid as the submitted autosomal sim:
#   0: p10m10  1: p10m30  2: p10m50    3: p20m10  4: p20m30  5: p20m50
#   6: p40m10  7: p40m30  8: p40m50    9: p60m10 10: p60m30 11: p60m50
#  12: p80m10 13: p80m30 14: p80m50
#
# NOTE: only e1, e2, e4 exist. e3 (in-trans amplification) and e5 (deletion on
# the other allele) are diploid-only constructs — see §20.1. Indices here are
# POSITIONAL, not the submitted exp numbers: exp_idx 2 is e4.

# --- machine-specific paths ---------------------------------------------------
# These used to be machine-specific literal defaults. On another machine those
# resolve to nothing, and a `${VAR:-/literal}` default means the script runs
# anyway against a path that does not exist. Paths now come from config.local.sh,
# which is per-machine and never committed. See README.md.
#
# The config is found by walking UP from this file. That is the one lookup that
# cannot itself be configured, so it must not depend on any absolute path.
# $VISOR_CONFIG overrides it, which is how a SLURM task that has been copied to a
# node-local spool directory passes the location through --export.
if [[ -z "${VISOR_CONFIG_LOADED:-}" ]]; then
    _cfg="${VISOR_CONFIG:-}"
    if [[ -z "$_cfg" ]]; then
        _d="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        while [[ "$_d" != "/" ]]; do
            [[ -r "$_d/config.local.sh" ]] && { _cfg="$_d/config.local.sh"; break; }
            _d="$(dirname "$_d")"
        done
    fi
    [[ -n "$_cfg" && -r "$_cfg" ]] || {
        echo "ERROR: config.local.sh not found by walking up from ${BASH_SOURCE[0]}." >&2
        echo "       cp config.example.sh config.local.sh at the repo root and edit it." >&2
        echo "       See README.md." >&2
        return 1 2>/dev/null || exit 1
    }
    # shellcheck disable=SC1090
    source "$_cfg"
    unset _cfg _d
fi

BASE_DIR="${BASE_DIR:-$CHRX_DIR}"
REF="${REF:-$REF_CHRX}"
AUTOSOME="${AUTOSOME:-chr22}"
HEMI="${HEMI:-chrX}"

# Reference pipeline scripts (externally owned and read-only here)
REPL_DIR="${REPL_DIR:-$REPLICATES_DIR}"
SEEDED_PY="${SEEDED_PY:-${REPL_DIR}/scripts/pipeline/visor_seeded.py}"

clone_genomes="${BASE_DIR}/truth/fastas/clone_genomes"

# --- Replicates ---------------------------------------------------------------
# REP selects an independent draw of SHORtS read noise, which is the quantity a
# block bootstrap has to resample: the submitted autosomal figure resamples its 30
# replicates for exactly this reason (the autosomal benchmark), and the
# chrX simulation had none.
#
# REP=0 is the ORIGINAL run and reproduces the previous seeds exactly. That is
# deliberate: the c25 results and their committed cn_bar/scores must stay
# reproducible, so adding replicate support cannot silently renumber them.
#
# Every replicate gets its own output tree as well as its own seed. A distinct
# seed with a shared path would have each replicate overwrite the last and leave
# one BAM set wearing the name of N runs -- the same failure the matched normal
# had when it was not coverage-named.
REP="${REP:-0}"
[[ "$REP" =~ ^[0-9]+$ ]] || { echo "ERROR: REP must be a non-negative integer, got '$REP'" >&2
                              return 1 2>/dev/null || exit 1; }
rep_tag=""; (( REP > 0 )) && rep_tag="_rep${REP}"

short_base="${BASE_DIR}/short${rep_tag}"

# config.local.sh exports SHORT_DIR and CALLS_DIR as FIXED paths, so they are
# overridden here rather than left to disagree with short_base: 11 reads
# $SHORT_DIR while 06 uses $short_base, and under REP=2 that combination would
# have genotyped replicate 0's BAMs and written them under replicate 2's name.
# Both are identical to the config values at REP=0.
SHORT_DIR="$short_base"
CALLS_DIR="${BASE_DIR}/calls${rep_tag}"

# Every other downstream directory is reached through a SEG_DIR / OUT_DIR
# override, so rather than enumerate them the consumers append $rep_tag to
# whatever they are given. SEG_DIR=$V/seg10k under REP=2 becomes
# $V/seg10k_rep2, which keeps the overrides rep-safe without a second knob.

# --- Coverage -----------------------------------------------------------------
# §21: SHORtS splits coverage by FASTA count, and chrX is present in only one of
# the two haplotypes, so chrX lands at HALF this value. COV=50 therefore gives
# chr22 50x and chrX 25x.
#
# COV IS REQUIRED, with no default: every path here is built from it, so a missing
# value is refused here rather than failing later on a missing file.
#
# 01_chrx_hack.sh builds haplotype FASTAs and has no coverage dimension at all,
# so it opts out rather than being made to pass a coverage it does not use.
if [[ -n "${CHRX_COV_OPTIONAL:-}" ]]; then
    cov="${COV:-}"
else
    [[ -n "${COV:-}" ]] || {
        echo "ERROR: COV is not set, and there is no default." >&2
        echo "       Every path below is derived from it. Run e.g. COV=50 $0" >&2
        return 1 2>/dev/null || exit 1
    }
    cov="$COV"
fi

# --- Aggregate outputs, coverage-named ----------------------------------------
# calls/, seg*/ and cn_bar*/ are keyed on <exp>_<cond>, and cond_name carries the
# coverage, so those are coverage-safe by construction. These three are NOT: they
# are one file per run, and they were fixed names. Regenerating the truth at a
# second coverage would therefore have overwritten the committed c25 truth table
# in place -- and 4,935 rows differing only in a label is not a difference anyone
# spots in a diff, while every c25 score already computed rested on it.
#
# Same failure as the coverage-blind matched normal, and it gets the same fix.
# Defined here rather than in each script so the readers cannot drift from the
# writer: 04 writes them, 09/12 read the conditions, 09/10/13 read the truth.
# Empty only under CHRX_COV_OPTIONAL, where they are not used either.
if [[ -n "$cov" ]]; then
    cond_tsv="${BASE_DIR}/truth/chrx_conditions_c${cov}.tsv"
    truth_tsv="${BASE_DIR}/truth/chrx_sv_ground_truth_c${cov}.tsv"
fi

# --- Bad nodes ----------------------------------------------------------------
# Nodes that kill tasks at start-up: exit 0:53 within seconds, and no log files
# written at all, which is the signature of a node that cannot reach /vast (SLURM
# fails to open the job's output files and gives up). It is not a pipeline error
# and it does not depend on the stage -- simulation, calling and segmentation
# have all been lost this way.
#
# Every node listed here has ZERO completed tasks against 1-12 failures in this
# run, so excluding them costs no capacity: they have never produced anything.
# Nodes with a mix of failures and completions are deliberately NOT listed --
# a task can fail for its own reasons, and dropping a working node would.
#
# Evidence, failures/completions over this run:
#   c592 80/0   c369 12/0   c632 8/0   c589 8/0   c586 4/0
#   c583  4/0   c581  4/0   c580 4/0   c577 4/0   c298 4/0
#
# Nodes go bad one at a time and mid-run: c592 was healthy when this list was
# first written and later took out 80 tasks across five replicates in one burst.
# So REGENERATE rather than trusting the literal below -- every node with
# failures and no completions, straight from the accounting:
#
#   sacct -S <start> -X -n -o State,NodeList,JobName%14 \
#     | awk '$3 ~ /chrx_/ && $2 ~ /^c[0-9]+$/ {
#         if ($1 ~ /FAILED|NODE_FAIL/) f[$2]++; else if ($1=="COMPLETED") c[$2]++ }
#       END { for (n in f) if (c[n]+0==0) printf "%s,", n }'
#
# c061 was on this list and has been REMOVED: over the full run it is 0 failures
# and 23 completions. It was added on a single early NODE_FAIL, which is the
# mistake this rule is meant to prevent -- a task can fail for its own reasons,
# and excluding a working node costs throughput for nothing. Require zero
# completions, not merely some failures.
#
# Override per invocation with CHRX_EXCLUDE_NODES=... , or set it empty to
# disable. This is a statement about one cluster on one day, not a permanent
# property of those machines.
CHRX_EXCLUDE_NODES="${CHRX_EXCLUDE_NODES:-}"

experiments=(e1 e2 e4)
purities=(10 10 10 20 20 20 40 40 40 60 60 60 80 80 80)
mixtures=(10 30 50 10 30 50 10 30 50 10 30 50 10 30 50)

decode_task() {
    local tid="$1"
    exp_idx=$(( tid / 15 ))
    cond_idx=$(( tid % 15 ))
    exp_name=${experiments[$exp_idx]}
    ppur=${purities[$cond_idx]}
    sub_mix=${mixtures[$cond_idx]}

    # Clone fractions: sub_mix% of tumour cells are c2, the rest c3, then normal
    # contamination. SVs shared by c2 and c3 are clonal; private SVs subclonal.
    c2_por=$(( ppur * sub_mix / 100 ))
    c3_por=$(( ppur - c2_por ))
    norm_por=$(( 100 - ppur ))
    frac=($c2_por $c3_por $norm_por)

    cond_name=c${cov}p${ppur}m${sub_mix}
    work_dir="${short_base}/${exp_name}/${cond_name}"

    # Unique per (replicate, experiment, condition). The 10000s digit was reserved
    # for the replicate index from the start; this spends it.
    #
    # REP=0 leaves the value identical to the pre-replicate formula, so the c25
    # run and the first c50 run stay reproducible from their recorded seeds.
    # Distinctness: exp_idx <= 2 and cond_idx <= 14 give a base of at most 3150,
    # well under the 10000 stride, so no two (REP, exp, cond) triples collide.
    visor_seed=$(( REP * 10000 + (exp_idx + 1) * 1000 + (cond_idx + 1) * 10 ))
}

# Build the SHORtS BED: max coordinate per contig across every source haplotype
# plus the reference. A haplotype shorter than the stated end is fine — pyfaidx
# truncates the slice and SHORtS scales read count to the sequence it actually
# retrieved, so per-base depth stays correct.
build_short_bed() {
    local out="$1"; shift
    cut -f1,2 "$@" "${REF}.fai" \
        | sort \
        | awk '$2 > maxvals[$1] {lines[$1]=$0; maxvals[$1]=$2} END {for (t in lines) print lines[t]}' \
        | awk 'OFS=FS="\t" {print $1, "1", $2, 100, 100}' \
        | sort -k1,1 > "$out"
}

# §21 GUARD. A contig absent from a haplotype is a WARNING with exit 0, and
# 00_visor_shorts.sh runs `set +o pipefail` precisely to tolerate VISOR warnings.
# So a contig-name typo would silently halve that contig's coverage and the run
# would look clean. Assert the warning count is exactly what we intend: one per
# clone directory, because every h2.fa deliberately lacks chrX.
check_warnings() {
    local log="$1" expected="$2"
    local seen
    seen=$(grep -c "Chromosome ${HEMI} not found" "$log" || true)
    local other
    other=$(grep -i 'warning' "$log" | grep -vc "Chromosome ${HEMI} not found" || true)
    printf '  %s-absent-from-h2 warnings: %d (expected %d)\n' "$HEMI" "$seen" "$expected"
    [[ "$seen" -eq "$expected" ]] || { echo "  GUARD FAILED: unexpected warning count" >&2; return 1; }
    [[ "$other" -eq 0 ]] || { echo "  GUARD FAILED: $other unexpected warning(s)" >&2
                              grep -i 'warning' "$log" | grep -v "Chromosome ${HEMI} not found" >&2; return 1; }
    return 0
}

# Mean depth at MAPQ>=20 over FIXED WINDOWS.
#
# WHOLE-CONTIG MEANS ARE NOT USABLE HERE (§29). `samtools depth -a` scores every
# reference position, including assembly N, at depth 0 — but SHORtS simulates
# reads only for non-N bases (`Nreads` uses `len(seq_) - Ns`). chr22 is 22.9% N
# (acrocentric p-arm) against chrX's 0.7%, so a whole-contig ratio reads
# (0.5 x 0.993)/0.771 ~= 0.64 for a perfectly correct hemizygous BAM. The
# original version of this guard would have failed a good normal.
#
# These windows are verified N-free. The chr22 four are the same AUTO_REF set
# the depth-based cn_bar estimators already use; the chrX four are
# clear of the exclusion bed (PAR, centromere, telomere, assembly gaps).
# 8 Mb per contig, equal on both sides.
AUTO_WIN_BED="${AUTO_WIN_BED:-${BASE_DIR}/resources/beds/depth_windows_auto.bed}"
HEMI_WIN_BED="${HEMI_WIN_BED:-${BASE_DIR}/resources/beds/depth_windows_hemi.bed}"

mean_depth() {
    local bam="$1" bed="$2"
    samtools depth -a -Q 20 -b "$bed" "$bam" \
        | awk '{s+=$3; n++} END {if (n>0) printf "%.4f\n", s/n; else print "NA"}'
}

depth_ratio() {
    local bam="$1" a x
    a=$(mean_depth "$bam" "$AUTO_WIN_BED")
    x=$(mean_depth "$bam" "$HEMI_WIN_BED")
    awk -v a="$a" -v x="$x" 'BEGIN {
        if (a == "NA" || x == "NA" || a + 0 == 0) { print "NA\tNA\tNA" }
        else { printf "%.4f\t%.3f\t%.3f\n", x / a, a, x } }'
}
