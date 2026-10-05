#!/bin/bash
#SBATCH --job-name=snp_vis
#SBATCH --output=log/02_snp_%A_%a.out
#SBATCH --error=log/02_snp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-2249%10
#
# GATK HaplotypeCaller + het SNP extraction near SVs.
# No scatter needed — VISOR uses only chr1-2 (~490Mb).

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
log_repro_header "02_snp_pipeline.sh"

mkdir -p "$snp_dir"

source "$CONDA_SH"
conda activate "$ENV_GATK"
echo "  gatk    : $(gatk --version 2>&1 | grep -v '^$' | head -1)"

# ---- GATK HaplotypeCaller (whole genome, chr1-2 only) ----
gatk --java-options "-Xmx6g" HaplotypeCaller \
    -R $ref \
    -I $tbam \
    -O $snp_dir/SNP.vcf.gz

bcftools view -v snps -g het -Oz -o $snp_dir/het_snp.vcf.gz $snp_dir/SNP.vcf.gz
tabix -p vcf $snp_dir/het_snp.vcf.gz

conda deactivate

# ---- SNP phasing on SV reads ----
conda activate "$ENV_VISOR"  # provides bcftools + samtools
echo "  bcftools: $(bcftools --version | head -1)"

# Get SV-flanking regions
Rscript "$get_sv_range" \
    -s "$cond_name" -o "$snp_dir" \
    -v "$svtyp_dir/svt_${cond_name}.vcf"

# Het SNPs near SVs
bcftools view -R $snp_dir/${cond_name}.bed \
    $snp_dir/het_snp.vcf.gz \
    -Ov -o $snp_dir/het_near_sv_${cond_name}.vcf

# Allele counts on SV-supporting reads
awk 'BEGIN{OFS="\t"} !/^#/ { print $1, $2-1, $2 }' \
    $snp_dir/het_near_sv_${cond_name}.vcf \
    > $snp_dir/pos_${cond_name}.bed

$samtools view -f 1 -F 2 -b $tbam > $snp_dir/sup_${cond_name}.bam
$samtools index $snp_dir/sup_${cond_name}.bam

bcftools mpileup \
    -f $ref -a DP,AD -A \
    -R $snp_dir/pos_${cond_name}.bed \
    $snp_dir/sup_${cond_name}.bam \
    -Ov > $snp_dir/het_on_sv_${cond_name}.vcf

rm $snp_dir/sup_${cond_name}.bam $snp_dir/sup_${cond_name}.bam.bai

conda deactivate
echo "Done: SNP pipeline [$rep / $exp_name / $cond_name]"
