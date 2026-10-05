#!/bin/bash
#SBATCH --job-name=pm_chrx_seg
#SBATCH --output=log/04a_chrx_seg_%A_%a.out
#SBATCH --error=log/04a_chrx_seg_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --array=0-329%50
#
# Stage 4a: chrX copy-number segmentation from read depth, for the prostate mixtures.
#
# WHY THIS STAGE EXISTS. The mixture subjects are male, so their X is single-copy and carries no
# heterozygous germline SNPs. FACETS derives copy number from those SNPs, so on X it returns "norm"
# for every segment -- the ABSENCE of a call, not evidence of copy-neutrality. Stage 04 therefore
# cannot supply the copy number SVCFit's hemizygous forms need, and this stage supplies it from
# depth instead.
#
# It drives helper/chrx_depth_segmentation.R, the depth-based cn_bar estimator.
#
# THE MIXTURES ARE GRCh37/hs37d5, NOT hg38. Their X contig is named "X" and is 155,270,560 bp; the
# cohort's is "chrX" at 156,040,895. The three CHRX_* variables below carry that difference. Left at
# the hg38 defaults, `samtools depth -r chrX:...` on a GRCh37 BAM returns no rows and the script
# would segment an empty chromosome rather than fail -- which is why the helper also refuses when every
# bin reads zero in both BAMs.
#
# psi_sample IS THE TUMOUR'S AUTOSOMAL PLOIDY, and it is not 2 here. FACETS reports 2.75 for 3m19
# rep1, so the psi_sample/2 factor that recovers absolute copy number matters far more than in the
# chrX simulation, where the value is 2 exactly by construction. It is read from the FACETS fit for
# this sample, so stage 04 must have run first.
#
# 11 conditions x 30 replicates = 330 tasks, indexed exactly as pipeline_common.sh does.

set -euo pipefail

_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"

# 4m AND 5m HAVE NO X. Visible in the data: their BAMs
# carry normal autosomal depth (20x, 30x) and EXACTLY ZERO on X, because 00_make_mixtures.sh builds
# the 4- and 5-clone mixtures from chromosome-split sources -- fil_bM_odd / fil_bM_even /
# fil_gM_odd / fil_gM_even (00_make_mixtures.sh:75-116). Those splits are autosomes only, so X
# belongs to neither pseudo-clone and drops out.
#
# Left to run, they do not fail cleanly: every bin reads zero, cn_bar is zero, log2(0) is -Inf, and
# DNAcopy dies inside smooth.CNA with "only 0's may be mixed with negative subscripts" -- an error
# that says nothing about the cause. Skipped here with exit 0, because an excluded condition is not
# a failed one and 60 red tasks in every future array would train the eye to ignore them.
#
# CONSEQUENCE FOR THE ANALYSIS: chrX rests on the 9 three-clone conditions, 9 x 30 = 270 samples,
# not 11 x 30 = 330. Say so wherever the chrX mixture results are reported.
case "$samp_name" in
    4m|5m)
        echo "SKIP: ${rep}/${samp_name} has no X by construction (chromosome-split mixture); nothing to segment"
        exit 0
        ;;
esac

# Repo-internal, so derived from the checkout rather than named absolutely.
# $_pm_dir is prostate_mixture/, set by pipeline_common.sh.
SEG_R="${CHRX_SEG_R:-$_pm_dir/scripts/helper/chrx_depth_segmentation.R}"
[[ -r "$SEG_R" ]] || { echo "ERROR: segmentation script not found: $SEG_R" >&2; exit 1; }

# DNAcopy is in the `visor` env and nowhere else on this machine -- not in svcfit, not in facet,
# which is where you would look first given this runs after stage 04. Resolved explicitly rather
# than left to whatever Rscript is on PATH, and checked, because the failure is a bare
# "no package called DNAcopy" three levels down in a 330-task array.
RSCRIPT="${CHRX_RSCRIPT:-${RSCRIPT}}"
[[ -x "$RSCRIPT" ]] || { echo "ERROR: Rscript not executable at $RSCRIPT" >&2; exit 1; }
"$RSCRIPT" -e 'if (!requireNamespace("DNAcopy", quietly=TRUE)) quit(status=1)' 2>/dev/null || {
    echo "ERROR: $RSCRIPT cannot load DNAcopy. Set CHRX_RSCRIPT to an R that has it." >&2; exit 1; }

# samtools is called by the segmentation script through system(); same reasoning.
export PATH="$(dirname "$SAMTOOLS"):$PATH"
[[ -x "$SAMTOOLS" ]] || { echo "ERROR: project samtools is not executable" >&2; exit 1; }

# GRCh37/hs37d5. Autosomal baseline regions are the cohort's, with the prefix dropped -- same loci,
# and GRCh37 coordinates differ from hg38 but these are 2 Mb windows in gene-poor regions chosen for
# a depth baseline, not for any feature, so the shift is immaterial.
export CHRX_CONTIG="${CHRX_CONTIG:-X}"
export CHRX_LEN="${CHRX_LEN:-155270560}"
export CHRX_AUTO_REF="${CHRX_AUTO_REF:-1:50000000-52000000,2:100000000-102000000,7:20000000-22000000,12:60000000-62000000}"
BIN_KB="${BIN_KB:-100}"

# WHICH TREE TO RUN AGAINST. pipeline_common.sh points at bootstrap/, the redo's output, which
# stage 00 has to create first. The PREVIOUS run is complete and on disk -- 30/30 replicates with
# all 11 mixture BAMs under replicates/ and all 11 FACETS fits under prostate_replicates/ -- so
# this stage can run there today without waiting.
#
# Set PM_BAM_ROOT and PM_WORK_ROOT together to point at it:
#   PM_BAM_ROOT="$MIX_BAM_DIR"
#   PM_WORK_ROOT="$WORK_ROOT"
# Unset, the pipeline_common.sh defaults apply and this runs against bootstrap/.
if [[ -n "${PM_BAM_ROOT:-}" ]];  then tbam="${PM_BAM_ROOT}/${rep}/${samp_name}.bam"; fi
if [[ -n "${PM_WORK_ROOT:-}" ]]; then
    work_dir="${PM_WORK_ROOT}/${rep}"
    fac_dir="${work_dir}/facet/${samp_name}"
fi

seg_dir="${work_dir}/chrx_seg"
out="${seg_dir}/chrx_seg_${samp_name}.csv"
mkdir -p "$seg_dir"

# psi_sample from this sample's FACETS fit. fac_dir ALREADY ends in $samp_name
# (pipeline_common.sh:30), so the file is ${fac_dir}/${samp_name}.RData -- not a further
# ${samp_name}/ below it.
facet_rdata="${fac_dir}/${samp_name}.RData"
[[ -s "$facet_rdata" ]] || { echo "ERROR: no FACETS fit at $facet_rdata -- run stage 04 first" >&2; exit 1; }
psi=$("$RSCRIPT" -e "e<-new.env(); load('${facet_rdata}', envir=e); f<-get('fit',envir=e); cat(f\$ploidy)" 2>/dev/null)
[[ -n "$psi" ]] || { echo "ERROR: could not read fit\$ploidy from $facet_rdata" >&2; exit 1; }

echo "========================================================"
echo "  task    : ${SLURM_ARRAY_TASK_ID:-NA}   ${rep} / ${samp_name}"
echo "  contig  : ${CHRX_CONTIG}   length: ${CHRX_LEN}   bin: ${BIN_KB} kb"
echo "  tumour  : ${tbam}"
echo "  normal  : ${nbam}"
echo "  psi     : ${psi}   (FACETS ploidy; NOT assumed to be 2)"
echo "  out     : ${out}"
echo "========================================================"

for f in "$tbam" "$nbam"; do
    [[ -s "$f" ]] || { echo "ERROR: missing BAM: $f" >&2; exit 1; }
done

"$RSCRIPT" "$SEG_R" "$tbam" "$nbam" "$psi" "$out" "$BIN_KB"

[[ -s "$out" ]] || { echo "ERROR: no segmentation written to $out" >&2; exit 1; }
echo "RESULT: ${rep}/${samp_name} -> ${out} ($(( $(wc -l < "$out") - 1 )) segments)"
