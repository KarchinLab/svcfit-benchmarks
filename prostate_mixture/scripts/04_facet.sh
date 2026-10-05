#!/bin/bash
#SBATCH --job-name=prst_facet
#SBATCH --output=log/04_facet_%A_%a.out
#SBATCH --error=log/04_facet_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --array=0-329%10
#
# Independent of 01-03; runs in parallel with SNP steps.

set -euo pipefail

# Locate pipeline_common.sh beside this script (same lookup as 04a/04b).
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"
mkdir -p $fac_dir $work_dir/facet

if [[ -f "$work_dir/facet/${samp_name}.bed" ]]; then
    echo "[$rep / $samp_name] Output already exists, skipping: ${samp_name}.bed"
    exit 0
fi

log_repro_header "04_facet.sh"
echo "  snp_vcf : ${snp_vcf}"
echo "  snp-pileup params: -g -q15 -Q20 -P100 -r25,0"

echo "[$rep / $samp_name] FACET"

: "${FACET_SNP_PILEUP:?set FACET_SNP_PILEUP in config.local.sh}"
: "${FACET_R_SCRIPT:?set FACET_R_SCRIPT in config.local.sh}"
[[ -x "$FACET_SNP_PILEUP" ]] || { echo "ERROR: FACET_SNP_PILEUP is not executable: $FACET_SNP_PILEUP" >&2; exit 1; }
[[ -f "$FACET_R_SCRIPT" ]] || { echo "ERROR: FACET_R_SCRIPT not found: $FACET_R_SCRIPT" >&2; exit 1; }

source "$PIPELINE_CONDA_SH"
conda activate "$ENV_FACET"
echo "  R       : $(Rscript --version 2>&1 | head -1)"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd $fac_dir
"$FACET_SNP_PILEUP" -g \
    -q15 -Q20 -P100 -r25,0 \
    $snp_vcf \
    ${samp_name}.csv \
    $nbam $tbam

Rscript "$FACET_R_SCRIPT" -s "$samp_name" -p "$fac_dir"

# Flat copy read by R analysis script
cp $fac_dir/${samp_name}.bed $work_dir/facet/${samp_name}.bed

conda deactivate
echo "Done: FACET [$rep / $samp_name]"
