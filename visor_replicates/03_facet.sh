#!/bin/bash
#SBATCH --job-name=facet_vis
#SBATCH --output=log/03_facet_%A_%a.out
#SBATCH --error=log/03_facet_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-2249%10
#
# FACET CNV segmentation on VISOR replicate BAMs.

set -euo pipefail

# pipeline_common.sh was sourced by absolute path, which only resolved on the
# machine it was written on. Locate it relative to this script instead;
# SLURM_SUBMIT_DIR first because sbatch copies the script to a node-local spool
# directory, where BASH_SOURCE points at the copy, not the checkout.
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"
log_repro_header "03_facet.sh"

mkdir -p $fac_dir

source "$CONDA_SH"
conda activate "$ENV_FACET"
echo "  R       : $(Rscript --version 2>&1 | head -1)"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd $fac_dir

# snp-pileup writes to tmp_${cond_name}.csv.gz; then strip 'chr' prefix for FACET
"$FACET_SNP_PILEUP" -g \
    -q15 -Q20 -P100 -r25,0 \
    $snp_dir/SNP.vcf.gz \
    tmp_${cond_name}.csv \
    $nbam $tbam

zcat tmp_${cond_name}.csv.gz \
    | sed '1!s/^chr//' \
    | gzip -c > ${cond_name}.csv.gz

rm tmp_${cond_name}.csv.gz

Rscript "$FACET_R_SCRIPT" -s "$cond_name" -p "$fac_dir"

conda deactivate
echo "Done: FACET [$rep / $exp_name / $cond_name]"
