#!/bin/bash
# Sourced by scripts 01-05 to decode SLURM_ARRAY_TASK_ID and set shared paths.
# 30 replicates x 11 conditions = 330 array tasks (0-329)
# cond_idx = task % 11   (0=3m19 ... 8=3m91, 9=4m, 10=5m)
# rep_idx  = task / 11   (0=rep1 ... 29=rep30)

conditions=(3m19 3m28 3m37 3m46 3m55 3m64 3m73 3m82 3m91 4m 5m)
cond_idx=$(( SLURM_ARRAY_TASK_ID % 11 ))
rep_idx=$(( SLURM_ARRAY_TASK_ID / 11 ))
samp_name=${conditions[$cond_idx]}
rep=rep$(( rep_idx + 1 ))   # rep1 … rep30

# --- machine-specific paths ---------------------------------------------------
# Two kinds of path are resolved two different ways, deliberately:
#
#   inside the repo   derived from THIS FILE's location, so they follow the checkout
#                     and need no configuration at all (ref, the helper scripts)
#   outside the repo  from config.local.sh, because they are a collaborator's source
#                     data that no checkout can locate (BAMs, intervals)
#
# The config is found by walking UP from this file, the same lookup
# visor_chrX/scripts/chrx_common.sh uses. See README.md.
_pm_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -z "${VISOR_CONFIG_LOADED:-}" ]]; then
    _cfg="${VISOR_CONFIG:-}"
    if [[ -z "$_cfg" ]]; then
        _d="$_pm_dir"
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

: "${PROSTATE_DATA_DIR:?set PROSTATE_DATA_DIR in config.local.sh -- the mixture source data (og_bam, split_chrom, bootstrap, intervals)}"
[[ -d "$PROSTATE_DATA_DIR" ]] || {
    echo "ERROR: PROSTATE_DATA_DIR -> $PROSTATE_DATA_DIR is not a directory" >&2
    return 1 2>/dev/null || exit 1
}

# Source BAMs — length-filtered copies (L=100, required by SVclone).
# Originals are untouched; filtered versions produced by 00a_filter_source_bams.sh.
src_bams=${PROSTATE_DATA_DIR}/prostate_replicates/source_bams
bM_bam=${src_bams}/fil_bM.bam
gM_bam=${src_bams}/fil_gM.bam
nbam=${PROSTATE_DATA_DIR}/og_bam/aWB_recal_sorted.bam
split=${src_bams}   # filtered split-chrom BAMs co-located with full BAMs
og_bam=${PROSTATE_DATA_DIR}/og_bam
split_chrom=${PROSTATE_DATA_DIR}/split_chrom

# GRCh37 reference (hs37d5). Defaults to prostate_mixture/hs37d5.fa; set PROSTATE_REF to use
# another copy. BWA indices are not shipped (about 5 GB; they regenerate in under an hour).
ref="${PROSTATE_REF:-$_pm_dir/hs37d5.fa}"

if [[ ! -f "$ref" ]]; then
    echo "ERROR: reference not found at $ref" >&2
    echo "       Set PROSTATE_REF, or place hs37d5.fa at $_pm_dir/" >&2
    return 1 2>/dev/null || exit 1
fi

# No index check here: the FASTA may be indexed in place before alignment.

# Mixture BAM produced by 00_make_mixtures.sh.
#
# Two trees:
#
#   replicates/rep*/           the mixture BAMs that step 00 writes and 01 reads
#   prostate_replicates/rep*/  everything 01-05 produce -- SNP, manta, svtyp, facet,
#                              svclone, cn_bar, chrx_seg -- and what the SVCFit stage reads
mix_bam_dir=${MIX_BAM_DIR:-${PROSTATE_DATA_DIR}/replicates}
work_root=${WORK_ROOT:-${PROSTATE_DATA_DIR}/prostate_replicates}

tbam=${PM_TBAM:-${mix_bam_dir}/${rep}/${samp_name}.bam}

# Analysis output directories
work_dir=${work_root}/$rep
manta_dir=$work_dir/manta/$samp_name
svtyp_dir=$work_dir/svtyp/$samp_name
snp_dir=$work_dir/SNP/$samp_name
fac_dir=$work_dir/facet/$samp_name
svc_dir=$work_dir/svclone/$samp_name

# Shared resources (intervals, scripts, config)
# snp_vcf: merged GATK HaplotypeCaller output from step 03 (same-sample SNPs, used by FACET snp-pileup)
snp_vcf=${snp_dir}/SNP.vcf.gz
interval_dir=${PROSTATE_DATA_DIR}/intervals

# The helper scripts are IN THE REPOSITORY and are resolved from this file's
# location, not from PROSTATE_DATA_DIR. They used to come out of the collaborator's
# tree, which meant the pipeline in the repo ran helpers that were not.
pipeline_dir=$_pm_dir
get_sv_range=${pipeline_dir}/scripts/helper/get_sv_range.R
make_input=${pipeline_dir}/scripts/helper/make_input.R
svclone_cfg=${pipeline_dir}/scripts/helper/svclone_config.ini

samtools="${SAMTOOLS:?set SAMTOOLS in config.local.sh}"
ppur=0.475

# ---------------------------------------------------------------
# Call at the top of every pipeline script (after sourcing this
# file) to write a reproducibility header to the SLURM log.
# The log files already capture stdout/stderr, so this is enough
# to reconstruct exactly what ran on any given task.
# ---------------------------------------------------------------
log_repro_header() {
    local script="${1:-unknown}"
    echo ""
    echo "========================================================"
    echo "  SCRIPT  : ${script}"
    echo "  DATE    : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  HOST    : $(hostname -s)"
    echo "  SLURM_JOB_ID        : ${SLURM_JOB_ID:-N/A}"
    echo "  SLURM_ARRAY_JOB_ID  : ${SLURM_ARRAY_JOB_ID:-N/A}"
    echo "  SLURM_ARRAY_TASK_ID : ${SLURM_ARRAY_TASK_ID:-N/A}"
    echo "  rep     : ${rep:-N/A}    condition : ${samp_name:-N/A}"
    echo "  cond_idx: ${cond_idx:-N/A}    rep_idx   : ${rep_idx:-N/A}"
    echo "--------------------------------------------------------"
    echo "  bM_bam  : ${bM_bam}"
    echo "  gM_bam  : ${gM_bam}"
    echo "  nbam    : ${nbam}"
    echo "  tbam    : ${tbam}"
    echo "  ref     : ${ref}"
    echo "  work_dir: ${work_dir}"
    echo "========================================================"
    echo ""
}
