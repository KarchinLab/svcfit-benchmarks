#!/bin/bash
#SBATCH --job-name=svclone_vis
#SBATCH --output=log/04_svclone_%A_%a.out
#SBATCH --error=log/04_svclone_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --array=0-2249%100
#
# SVclone on VISOR replicate BAMs.

set -euo pipefail

# pipeline_common.sh was sourced by absolute path, which only resolved on the
# machine it was written on. Locate it relative to this script instead;
# SLURM_SUBMIT_DIR first because sbatch copies the script to a node-local spool
# directory, where BASH_SOURCE points at the copy, not the checkout.
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"
log_repro_header "04_svclone.sh"

sv_file=$svtyp_dir/svt_${cond_name}.vcf
ppur_frac=$(echo "scale=3; $ppur / 100" | bc)

# Guard rm -rf
if [[ -z "$svc_dir" || "$svc_dir" != "$work_dir"/* ]]; then
    echo "ERROR: svc_dir is empty or not under work_dir: '$svc_dir'" >&2; exit 1
fi
rm -rf "$svc_dir"
mkdir -p "$svc_dir"
cd "$svc_dir"

export PATH="$ENV_SVCLONE/bin:$PATH"
SVCLONE="$ENV_SVCLONE/bin/svclone"
SVCLONE_RSCRIPT="$ENV_SVCLONE/bin/Rscript"
[[ -x "$SVCLONE" && -x "$SVCLONE_RSCRIPT" ]] || {
    echo "ERROR: incomplete configured SVclone environment: $ENV_SVCLONE" >&2; exit 1;
}
echo "  svclone : $("$SVCLONE" --version 2>&1 | head -1)"
echo "  ppur    : ${ppur_frac}"
echo "  sv_file : ${sv_file}"

"$SVCLONE_RSCRIPT" "$make_input" \
    -p "$ppur_frac" -f "$sv_file" \
    -s "$cond_name" -o "$svc_dir" -c "$fac_dir"

echo 'Annotating...'
"$SVCLONE" annotate \
    -i svc_${cond_name}_simple.txt \
    -b "$tbam" -s "$cond_name" \
    --sv_format simple -cfg "$svclone_cfg"

echo 'Counting...'
"$SVCLONE" count \
    -i "${cond_name}/${cond_name}_svin.txt" \
    -b "$tbam" -s "$cond_name" -cfg "$svclone_cfg"

echo 'Filtering...'
"$SVCLONE" filter \
    -s "$cond_name" \
    -i "${cond_name}/${cond_name}_svinfo.txt" \
    -p "spp_${cond_name}.txt" \
    -c cnv.txt -cfg "$svclone_cfg"

echo 'Clustering...'
"$SVCLONE" cluster -s "$cond_name" -cfg "$svclone_cfg"

rdata="${cond_name}/ccube_out/${cond_name}_ccube_sv_results.RData"
if [[ -f "$rdata" ]]; then
    "$SVCLONE_RSCRIPT" -e 'load(commandArgs(TRUE)[1]); stopifnot(exists("doubleBreakPtsRes"))' "$rdata"
    echo "Done: SVclone [$rep / $exp_name / $cond_name]"
else
    echo "ERROR: $rdata not produced" >&2; exit 1
fi
