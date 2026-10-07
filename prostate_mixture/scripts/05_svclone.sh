#!/bin/bash
#SBATCH --job-name=prst_svclone
#SBATCH --output=log/05_svclone_%A_%a.out
#SBATCH --error=log/05_svclone_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --array=0-329%100
#
# Depends on 01 (Manta VCF) and 04 (FACET BED).

set -euo pipefail

# Locate pipeline_common.sh beside this script (same lookup as 04a/04b).
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"

# SV input is the SVtyper VCF from step 01 (same records as the Manta VCF, re-genotyped), the same file SVCFit reads.
sv_file=${PM_SV_FILE:-$svtyp_dir/svt_${samp_name}.vcf}
fac_input=${PM_FAC_DIR:-$fac_dir}
rdata=${svc_dir}/${samp_name}/ccube_out/${samp_name}_ccube_sv_results.RData

# 4m and 5m are built from chromosome-split, autosome-only sources: no chrX, so the chrX rerun has
# nothing to add. Skip them BEFORE the rm -rf so their existing output is preserved untouched.
if [[ "$samp_name" == "4m" || "$samp_name" == "5m" ]]; then
    echo "[$rep / $samp_name] no chrX by construction, skipping"
    exit 0
fi

# FORCE=1 regenerates even when the RData exists. The chrX rerun MUST force, because the existing
# RData is the autosome-only output and would otherwise short-circuit here.
if [[ -f "$rdata" && "${FORCE:-0}" != "1" ]]; then
    echo "[$rep / $samp_name] Output already exists, skipping: $rdata"
    exit 0
fi

# Guard: only rm -rf if svc_dir is a non-empty, absolute path under work_dir
if [[ -z "$svc_dir" || "$svc_dir" != "$work_dir"/* ]]; then
    echo "ERROR: svc_dir is empty or not under work_dir: '$svc_dir'" >&2; exit 1
fi
rm -rf "$svc_dir"
mkdir -p "$svc_dir"
cd $svc_dir

log_repro_header "05_svclone.sh"
echo "  ppur    : ${ppur}"
echo "  config  : ${svclone_cfg}"
echo "  sv_file : ${sv_file}"
echo "  fac_input: ${fac_input}"
echo "  make_input: ${make_input}"

echo "[$rep / $samp_name] SVclone"

# Put the configured project prefix first for both this script and child processes
# launched by SVclone.  Do not rely on `conda activate`: the core runtime is also
# configured explicitly on PATH, so activation cannot reliably reorder it.
# SVclone needs Python 3.10, numpy<2, setuptools<81, pandas<3, and R optparse/ccube;
# these are pinned and checked by tools/setup_rockfish_tool_envs.sh.
export PATH="$ENV_SVCLONE/bin:$PATH"
SVCLONE="$ENV_SVCLONE/bin/svclone"
SVCLONE_RSCRIPT="$ENV_SVCLONE/bin/Rscript"
echo "  svclone : $($SVCLONE --version 2>&1 | head -1)"

$SVCLONE_RSCRIPT $make_input \
    -p $ppur -f $sv_file \
    -s $samp_name -o $svc_dir -c $fac_input

echo 'Annotating...'
$SVCLONE annotate \
    -i svc_${samp_name}_simple.txt \
    -b $tbam -s $samp_name \
    --sv_format simple -cfg $svclone_cfg

echo 'Counting...'
$SVCLONE count \
    -i ${samp_name}/${samp_name}_svin.txt \
    -b $tbam -s $samp_name -cfg $svclone_cfg

echo 'Filtering...'
$SVCLONE filter \
    -s $samp_name \
    -i ${samp_name}/${samp_name}_svinfo.txt \
    -p spp_${samp_name}.txt \
    -c cnv.txt -cfg $svclone_cfg

echo 'Clustering...'
$SVCLONE cluster -s $samp_name -cfg $svclone_cfg

if [ -f $rdata ]; then
    echo "Done: SVclone [$rep / $samp_name]"
else
    echo "ERROR: $rdata not produced" >&2; exit 1
fi
