#!/bin/bash
#SBATCH --job-name=pm_chrx_join
#SBATCH --output=log/04b_chrx_join_%A_%a.out
#SBATCH --error=log/04b_chrx_join_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:30:00
#SBATCH --array=0-329%50
#
# Stage 4b: segment-to-SV join for the prostate mixtures. Turns 04a's per-segment cn_bar into the
# per-SV cn_bar (and bg_cn) that SVCFit's hemizygous forms consume.
#
# Drives visor_chrX/scripts/08_chrx_segment_sv_join.py -- the same program the simulation and the
# cohort use -- so all three share one join rule (maximum-overlap segment, taken wholesale; see that
# file's header).
#
# WHERE THE SV LIST COMES FROM, because the three pipelines differ and this one had no precedent:
#   simulation  truth/chrx_sv_ground_truth_c*.tsv -- truth is known there
#   cohort      a prior SVCFit run's bed, chrX rows
#   mixtures    the SVtyper VCF, here
# The mixtures run SVCFit inside the step-06 Rmd and produce no per-sample bed, so there is no bed
# to read. The called SVs are the right list regardless: these are the variants SVCFit will estimate.
#
# BREAKENDS HAVE NO END, and are not dropped. 4 of the 14 X records in rep1/3m19 are BND, where
# Manta emits no END= in INFO. A breakend is a single position, so start = end = POS and the join
# assigns it whichever segment contains that point. This is exactly what the cohort does
# (`e = $h["END"]; if (e=="NA" || e=="") e = $h["POS"]`), copied rather than reinvented so all three
# pipelines treat breakends identically.
#
# GRCh37/hs37d5: contig "X", length 155,270,560. The VCFs use bare names -- verified in both the
# Manta and SVtyper output, which also carry the hs37d5 decoy contig.

set -euo pipefail

_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"

# See 04a: these two have no X by construction and never will.
case "$samp_name" in
    4m|5m)
        echo "SKIP: ${rep}/${samp_name} has no X by construction (chromosome-split mixture)"
        exit 0
        ;;
esac

# Repo-internal; see the note in 04a. $_pm_dir is prostate_mixture/.
JOIN_PY="${JOIN_PY:-$_pm_dir/../visor_chrX/scripts/08_chrx_segment_sv_join.py}"
[[ -r "$JOIN_PY" ]] || { echo "ERROR: join script not found: $JOIN_PY" >&2; exit 1; }

CHROM="${CHRX_CONTIG:-X}"
CHROM_LEN="${CHRX_LEN:-155270560}"
# WRITE chrX, NOT X. SVCFit's load_data() rewrites every CHROM to
# paste0("chr", sub("^chr","",CHROM)), so although these VCFs say "X", calc_svcf matches
# hemi_cn_bar on (CHROM, POS) with CHROM = "chrX". A table saying "X" would miss every row and
# all 2714 chrX SVs would fall through as having no depth -- silently, since an unmatched
# lookup is indistinguishable from an SV that genuinely had no segment.
# The VCF is still FILTERED on "X" (--chrom); only the output name changes.
OUT_CHROM="${CHRX_OUT_CONTIG:-chrX}"

# Same tree selection as 04a.
if [[ -n "${PM_BAM_ROOT:-}" ]];  then tbam="${PM_BAM_ROOT}/${rep}/${samp_name}.bam"; fi
if [[ -n "${PM_WORK_ROOT:-}" ]]; then work_dir="${PM_WORK_ROOT}/${rep}"; fi

seg="${work_dir}/chrx_seg/chrx_seg_${samp_name}.csv"
vcf="${work_dir}/svtyp/${samp_name}/svt_${samp_name}.vcf"
OUT="${work_dir}/cn_bar"
mkdir -p "$OUT"

[[ -s "$seg" ]] || { echo "ERROR: no segmentation at $seg -- run 04a first" >&2; exit 1; }
[[ -s "$vcf" ]] || { echo "ERROR: no SVtyper VCF at $vcf -- run stage 01 first" >&2; exit 1; }

tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

# VCF -> CHROM/start/end/ID, X only. END from INFO, falling back to POS for breakends.
awk -F'\t' -v OFS='\t' -v C="$CHROM" '
    BEGIN { print "CHROM","start","end","ID" }
    /^#/ { next }
    $1 == C {
        end = ""
        if (match($8, /(^|;)END=[0-9]+/)) { s = substr($8, RSTART, RLENGTH); sub(/^;/, "", s); end = substr(s, 5) }
        if (end == "") end = $2
        print $1, $2, end, $3
    }' "$vcf" > "$tmp/svs.tsv"

n=$(( $(wc -l < "$tmp/svs.tsv") - 1 ))
nseg=$(( $(wc -l < "$seg") - 1 ))
echo "========================================================"
echo "  ${rep} / ${samp_name}   contig ${CHROM}   len ${CHROM_LEN}"
echo "  SVs on ${CHROM}: ${n}    segments: ${nseg}"
echo "========================================================"

# A sample with no X SVs is a real outcome, not an error -- but it must not produce an empty
# cn_bar table that reads downstream as "joined, nothing to resolve". Skip it explicitly.
if (( n == 0 )); then
    echo "NOTE: ${samp_name} has no SVs on ${CHROM}; nothing to join"
    exit 0
fi

python3 "$JOIN_PY" \
    --segments "$seg" --svs "$tmp/svs.tsv" \
    --chrom-col CHROM --start-col start --end-col end \
    --chrom "$CHROM" --chrom-len "$CHROM_LEN" --out-chrom "$OUT_CHROM" \
    --out    "${OUT}/${samp_name}_cn_bar.csv" \
    --bg-out "${OUT}/${samp_name}_bg_cn.csv" \
    --audit  "${OUT}/${samp_name}_audit.tsv" \
    ${MODE_FLAG:+$MODE_FLAG} --quiet

[[ -s "${OUT}/${samp_name}_cn_bar.csv" ]] || { echo "ERROR: no cn_bar table written" >&2; exit 1; }
[[ -s "${OUT}/${samp_name}_bg_cn.csv" ]] || { echo "ERROR: no bg_cn table written" >&2; exit 1; }

# THE POS ROUND-TRIP, copied from the cohort for the same reason. calc_svcf joins cn_bar on
# (CHROM, POS), so a table that has lost or altered a position does not fail there -- it leaves
# that SV without depth, and the row falls back to unresolved, which is indistinguishable from an
# SV that genuinely had no segment. Compare the sets both ways.
awk -F'\t' 'NR>1{print $2}' "$tmp/svs.tsv" | sort -u > "$tmp/want"
awk -F',' 'NR>1{gsub(/"/,"",$2); print $2}' "${OUT}/${samp_name}_cn_bar.csv" | sort -u > "$tmp/got"
if ! diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ERROR: ${samp_name}: cn_bar POS set does not match the VCF's ${CHROM} POS set" >&2
    echo "  in VCF not in cn_bar:" >&2; comm -23 "$tmp/want" "$tmp/got" | head -5 | sed 's/^/    /' >&2
    echo "  in cn_bar not in VCF:" >&2; comm -13 "$tmp/want" "$tmp/got" | head -5 | sed 's/^/    /' >&2
    exit 1
fi

echo "RESULT: ${rep}/${samp_name} -> ${OUT}/${samp_name}_cn_bar.csv ($(( $(wc -l < "${OUT}/${samp_name}_cn_bar.csv") - 1 )) rows), bg_cn written, POS round-trip clean"
