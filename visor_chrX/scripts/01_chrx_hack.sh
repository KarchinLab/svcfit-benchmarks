#!/bin/bash
#SBATCH --job-name=chrx_hack
#SBATCH --output=log/chrx_hack_%j.out
#SBATCH --error=log/chrx_hack_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#
# Build hemizygous chrX clone genomes for the SVCFit chrX simulation.
# Hemizygous analogue of scripts/genome_setup/run_m_hack.sh (session log §20, §21).
#
# KARYOTYPE PER CLONE (verified empirically in §21):
#   h1.fa = chr22 + chrX(variants)     h2.fa = chr22 only
# SHORtS splits coverage by FASTA count (SHORtS.py:1039-1041) and skips a contig
# absent from a haplotype with a warning (:244-247). So chr22 draws from both
# haplotypes and chrX from one -> chrX sits at exactly half chr22's depth.
# A SINGLE fasta holding both contigs would make chr22 haploid-equivalent and
# silently break psi_sample = 2. The two-FASTA split is required.
#
# EXPERIMENTS: e1, e2, e4 only. e3 (in-trans amp) and e5 (deletion on the other
# allele) are diploid-only constructs — see §20.1.
#
# NO SNP BACKGROUND anywhere: a male chrX carries no het
# germline SNPs, so clone genomes are built straight off plain reference. Spiking
# SNPs would let SVCFit form an allele copy ratio and validate a code path real
# data never takes, while looking normal.

set -euo pipefail

# BASE_DIR / REF / AUTOSOME / HEMI were literal defaults here, duplicating
# chrx_common.sh. Both copies had to be edited per machine, and only one ever was.
# chrx_common.sh now resolves them from config.local.sh; see README.md.
# This stage builds haplotype FASTAs and has no coverage dimension, so it opts
# out of the COV requirement rather than passing a coverage it never uses.
CHRX_COV_OPTIONAL=1
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/chrx_common.sh"

sv_beds="${BASE_DIR}/truth/sv_beds"
out="${BASE_DIR}/truth/fastas/clone_genomes"
work="${BASE_DIR}/truth/fastas/.work"

# NOTE ON LAYOUT: §12(b) recorded that run_m_hack.sh writes truth/ while
# pipeline_common.sh reads ground_truth/. This standalone chrX simulation uses
# truth/ throughout, matching publication_package/'s real layout.

command -v VISOR   >/dev/null || { echo "ERROR: VISOR not on PATH from ENV_VISOR" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }
[[ -s "${REF}.fai" ]] || { echo "ERROR: missing ${REF}.fai" >&2; exit 1; }

REF_AUTO_LEN=$(awk -v c="$AUTOSOME" '$1==c{print $2}' "${REF}.fai")
REF_HEMI_LEN=$(awk -v c="$HEMI"     '$1==c{print $2}' "${REF}.fai")
[[ -n "$REF_AUTO_LEN" && -n "$REF_HEMI_LEN" ]] || {
    echo "ERROR: ${REF}.fai lacks $AUTOSOME and/or $HEMI" >&2; exit 1; }

echo "=========================================================="
echo "  reference : $REF"
echo "  $AUTOSOME : $REF_AUTO_LEN   $HEMI : $REF_HEMI_LEN"
echo "  base      : $BASE_DIR"
echo "  host      : $(hostname -s)   job: ${SLURM_JOB_ID:-interactive}"
echo "=========================================================="

rm -rf "$work"; mkdir -p "$work" "$out"

# ---------------------------------------------------------------------------
# Expected length change of a contig after applying a HACk bed.
#   deletion            -L        inversion  0        tandem duplication +L*(copies-1)
# BED is half-open here, so L = end - start. Reported against the observed value
# with a tolerance of 1 bp per variant to absorb end-convention differences; a
# genuinely dropped variant is 100 kb off and will trip this.
# ---------------------------------------------------------------------------
expected_delta() {
    awk -F'\t' '
        BEGIN                      { d = 0; n = 0 }
        $4 == "deletion"           { d -= ($3 - $2) }
        $4 == "inversion"          { }
        $4 == "tandem duplication" { d += ($3 - $2) * ($5 - 1) }
        { d += $6; n++ }   # col 6 = breakseqlen: VISOR inserts this many random
                           # bases at each breakpoint. Worth 354 bp on c2.bed, so
                           # omitting it made this check false-alarm (§23).
        END { print d "\t" n }
    ' "$1"
}

# h2 is identical for every clone: the untouched autosome, no variants.
echo "[$(date +%T)] building shared h2 (${AUTOSOME} only)"
samtools faidx "$REF" "$AUTOSOME" > "$work/h2.fa"
samtools faidx "$work/h2.fa"

# ---------------------------------------------------------------------------
# assemble <exp> <clone> <final h1.fa produced by HACk>
# ---------------------------------------------------------------------------
assemble() {
    local exp="$1" cln="$2" h1src="$3"
    local dest="$out/$exp/$cln"
    mkdir -p "$dest"
    mv "$h1src" "$dest/h1.fa"
    cp "$work/h2.fa" "$dest/h2.fa"
    samtools faidx "$dest/h1.fa"
    samtools faidx "$dest/h2.fa"
}

# ---------------------------------------------------------------------------
# Clone genomes
#   c2 / c3   carry the exp1 SV set        (c2.bed,  c3.bed)
#   c22 / c33 carry the CNV-overlap SV set (c22.bed, c33.bed) used by e2 and e4
# ---------------------------------------------------------------------------
clones=(c2 c3)
ccs=(c22 c33)

for i in 0 1; do
    cln=${clones[$i]}
    cc=${ccs[$i]}

    # --- e1: SV alone, copy-neutral --------------------------------------
    echo "[$(date +%T)] e1/$cln  <- $cln.bed"
    rm -rf "$work/e1_$cln"
    VISOR HACk -g "$REF" -b "$sv_beds/$cln.bed" -o "$work/e1_$cln"
    assemble e1 "$cln" "$work/e1_$cln/h1.fa"

    # --- e2: SV FIRST, then in-cis amplification -------------------------
    # SVCF = cn_bar*VAF - (cn_bar-1): the SV predates the amplification, so the
    # amplification carries it and all copies bear the SV.
    # Uses dup_${cc}_postsv.bed, NOT dup_${cc}.bed: the amp lands on a genome the
    # SV pass already reshaped, so it must be addressed in post-SV coordinates.
    # Reference coordinates here misplace the amplification and, for c33, overrun
    # the shortened contig outright (§23).
    echo "[$(date +%T)] e2/$cln  <- $cc.bed then dup_${cc}_postsv.bed"
    rm -rf "$work/e2a_$cln" "$work/e2_$cln"
    VISOR HACk -g "$REF"                  -b "$sv_beds/$cc.bed"            -o "$work/e2a_$cln"
    VISOR HACk -g "$work/e2a_$cln/h1.fa"  -b "$sv_beds/dup_${cc}_postsv.bed" -o "$work/e2_$cln"
    assemble e2 "$cln" "$work/e2_$cln/h1.fa"

    # --- e4: amplification FIRST, then SV --------------------------------
    # SVCF = cn_bar*VAF: the SV lands in one copy of an already-amplified region.
    # aft_*.bed coordinates are pre-shifted for the amplification (§20.4).
    echo "[$(date +%T)] e4/$cln  <- dup_$cc.bed then aft_$cc.bed"
    rm -rf "$work/e4a_$cln" "$work/e4_$cln"
    VISOR HACk -g "$REF"                  -b "$sv_beds/dup_$cc.bed" -o "$work/e4a_$cln"
    VISOR HACk -g "$work/e4a_$cln/h1.fa"  -b "$sv_beds/aft_$cc.bed" -o "$work/e4_$cln"
    assemble e4 "$cln" "$work/e4_$cln/h1.fa"
done

# ---------------------------------------------------------------------------
# Germline / normal-contamination clone: plain reference, no variants.
# Used BOTH as the matched normal and as the normal-contamination clone inside
# the tumour BAM. §12(a): both are the same trap with the same invisible failure
# mode — a diploid normal against a haploid tumour halves every cn_bar and every
# chrX fraction comes out 2x wrong with nothing appearing broken.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] normal (no variants)"
mkdir -p "$out/normal"
cp "$REF" "$out/normal/h1.fa"
cp "$work/h2.fa" "$out/normal/h2.fa"
samtools faidx "$out/normal/h1.fa"
samtools faidx "$out/normal/h2.fa"

# ---------------------------------------------------------------------------
# GUARDS. §21: SHORtS reports a missing contig as a warning and exits 0, and
# 00_visor_shorts.sh runs `set +o pipefail` precisely to tolerate warnings — so
# a wrong contig set here would silently halve coverage downstream. Verify the
# karyotype structurally now rather than trusting exit status later.
# ---------------------------------------------------------------------------
echo
echo "=== verification ==="
fail=0
printf '%-12s %-26s %-26s %s\n' DIR H1_CONTIGS H2_CONTIGS CHECKS

for d in "$out"/e{1,2,4}/c{2,3} "$out/normal"; do
    [[ -d $d ]] || continue
    label="${d#$out/}"
    h1c=$(cut -f1 "$d/h1.fa.fai" | sort | tr '\n' ',')
    h2c=$(cut -f1 "$d/h2.fa.fai" | sort | tr '\n' ',')
    notes=()

    # h1 must carry both contigs; h2 must carry the autosome ALONE
    # (sort puts chr22 before chrX: '2' < 'X')
    [[ "$h1c" == "${AUTOSOME},${HEMI}," ]] || notes+=("h1-contigs!")
    [[ "$h2c" == "${AUTOSOME}," ]] || notes+=("h2-contigs!")

    # the autosome is the depth baseline: it must be byte-for-byte the same length
    a1=$(awk -v c="$AUTOSOME" '$1==c{print $2}' "$d/h1.fa.fai")
    a2=$(awk -v c="$AUTOSOME" '$1==c{print $2}' "$d/h2.fa.fai")
    [[ "$a1" == "$REF_AUTO_LEN" ]] || notes+=("h1-$AUTOSOME-len=$a1!")
    [[ "$a2" == "$REF_AUTO_LEN" ]] || notes+=("h2-$AUTOSOME-len=$a2!")

    x1=$(awk -v c="$HEMI" '$1==c{print $2}' "$d/h1.fa.fai")
    [[ ${#notes[@]} -eq 0 ]] && notes=("ok")
    [[ "${notes[0]}" == "ok" ]] || fail=1
    printf '%-12s %-26s %-26s %s  (%s=%s)\n' "$label" "$h1c" "$h2c" "${notes[*]}" "$HEMI" "$x1"
done

echo
echo "=== ${HEMI} length: observed vs predicted from the beds ==="
check_len() {
    local label="$1" dir="$2"; shift 2
    local obs delta=0 nvar=0 r
    obs=$(awk -v c="$HEMI" '$1==c{print $2}' "$out/$dir/h1.fa.fai")
    for b in "$@"; do
        r=$(expected_delta "$sv_beds/$b")
        delta=$(( delta + $(echo "$r" | cut -f1) ))
        nvar=$(( nvar + $(echo "$r" | cut -f2) ))
    done
    local exp_len=$(( REF_HEMI_LEN + delta )) diff
    diff=$(( obs - exp_len )); diff=${diff#-}
    if [[ $diff -le $nvar ]]; then
        printf '  %-10s observed %-12d predicted %-12d diff %-6d OK (tol %d)\n' \
               "$label" "$obs" "$exp_len" "$diff" "$nvar"
    else
        printf '  %-10s observed %-12d predicted %-12d diff %-6d MISMATCH (tol %d)\n' \
               "$label" "$obs" "$exp_len" "$diff" "$nvar"
        fail=1
    fi
}

check_len e1/c2 e1/c2 c2.bed
check_len e1/c3 e1/c3 c3.bed
check_len e2/c2 e2/c2 c22.bed dup_c22_postsv.bed
check_len e2/c3 e2/c3 c33.bed dup_c33_postsv.bed
check_len e4/c2 e4/c2 dup_c22.bed aft_c22.bed
check_len e4/c3 e4/c3 dup_c33.bed aft_c33.bed

nx=$(awk -v c="$HEMI" '$1==c{print $2}' "$out/normal/h1.fa.fai")
printf '  %-10s observed %-12d predicted %-12d diff %-6d %s\n' \
       normal "$nx" "$REF_HEMI_LEN" "$(( nx - REF_HEMI_LEN ))" \
       "$([[ $nx -eq $REF_HEMI_LEN ]] && echo OK || { echo MISMATCH; fail=1; })"

echo
if [[ $fail -ne 0 ]]; then
    echo "RESULT: FAILED — do not proceed to SHORtS" >&2
    exit 1
fi
rm -rf "$work"
echo "RESULT: all clone genomes built and verified"
echo "next: SHORtS over $out/{e1,e2,e4}/{c2,c3} + normal, e.g. COV=50 (§21)"
