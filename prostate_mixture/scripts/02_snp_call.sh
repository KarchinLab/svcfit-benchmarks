#!/bin/bash
#SBATCH --job-name=prst_gatk
#SBATCH --output=log/02_snp_%A_%a.out
#SBATCH --error=log/02_snp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=40:00:00
#SBATCH --array=0-6929%30
#
# 6930 tasks = 330 (rep x condition) x 21 (GATK scatter intervals, 0000-0020)

set -euo pipefail

orig_task_id=$SLURM_ARRAY_TASK_ID
outer=$(( SLURM_ARRAY_TASK_ID / 21 ))
interval_idx=$(( SLURM_ARRAY_TASK_ID % 21 ))

export SLURM_ARRAY_TASK_ID=$outer
# Locate pipeline_common.sh beside this script (same lookup as 04a/04b).
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"

mkdir -p $snp_dir
mapfile -t interval_lst < <(find "$interval_dir" -name '*.interval_list' | sort)
interval=${interval_lst[$interval_idx]}

if [[ -f "$snp_dir/int_${interval_idx}_het_snp.vcf.gz.tbi" ]]; then
    echo "[$rep / $samp_name] Interval $interval_idx already done, skipping."
    exit 0
fi

log_repro_header "02_snp_call.sh"
echo "  orig_task_id: ${orig_task_id}  (outer=${outer}, interval_idx=${interval_idx})"
echo "  interval file: ${interval}"

source "$PIPELINE_CONDA_SH"
conda activate "$ENV_GATK"
echo "  gatk    : $(gatk --version 2>&1 | grep -v '^$' | head -1)"

echo "[$rep / $samp_name] GATK interval $interval_idx"

gatk --java-options "-Xmx5g" HaplotypeCaller \
    -R $ref \
    -I $tbam \
    -O $snp_dir/int_${interval_idx}_SNP.vcf.gz \
    -L $interval

bcftools view -v snps -g het -Oz \
    -o $snp_dir/int_${interval_idx}_het_snp.vcf.gz \
    $snp_dir/int_${interval_idx}_SNP.vcf.gz
tabix -p vcf $snp_dir/int_${interval_idx}_het_snp.vcf.gz

conda deactivate
echo "Done: GATK interval $interval_idx [$rep / $samp_name]"
