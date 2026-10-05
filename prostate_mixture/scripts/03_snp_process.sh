#!/bin/bash
#SBATCH --job-name=prst_snpproc
#SBATCH --output=log/03_snpproc_%A_%a.out
#SBATCH --error=log/03_snpproc_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=4:00:00
#SBATCH --array=0-329%10
#
# Depends on 01 (SVtyper VCF) and 02 (GATK intervals).

set -euo pipefail

# Locate pipeline_common.sh beside this script (same lookup as 04a/04b).
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"
mkdir -p $work_dir/SNP   # flat dir read by R script

if [[ -f "$work_dir/SNP/het_on_sv_${samp_name}.vcf" ]]; then
    echo "[$rep / $samp_name] Output already exists, skipping: het_on_sv_${samp_name}.vcf"
    exit 0
fi

log_repro_header "03_snp_process.sh"

echo "[$rep / $samp_name] SNP merge + extraction"

source "$PIPELINE_CONDA_SH"
conda activate "$ENV_GATK"
echo "  gatk    : $(gatk --version 2>&1 | grep -v '^$' | head -1)"

# Merge 30 GATK interval VCFs
gatk GatherVcfs \
    $(for i in $(seq 0 20); do echo "-I $snp_dir/int_${i}_SNP.vcf.gz"; done) \
    -O $snp_dir/SNP.vcf.gz
tabix $snp_dir/SNP.vcf.gz

gatk GatherVcfs \
    $(for i in $(seq 0 20); do echo "-I $snp_dir/int_${i}_het_snp.vcf.gz"; done) \
    -O $snp_dir/het_SNP_merged.vcf.gz
tabix $snp_dir/het_SNP_merged.vcf.gz

conda deactivate
source "$PIPELINE_CONDA_SH"
conda activate "$ENV_VISOR"
echo "  bcftools: $(bcftools --version | head -1)"
echo "  mpileup params: -f $ref -a DP -A"

# SV-flanking windows from SVtyper VCF
Rscript $get_sv_range \
    -s $samp_name -o $snp_dir \
    -v $svtyp_dir/svt_${samp_name}.vcf

# Het SNPs near SVs
bcftools view -R $snp_dir/${samp_name}.bed \
    $snp_dir/het_SNP_merged.vcf.gz \
    -Ov -o $snp_dir/het_near_sv_${samp_name}.vcf

# Discordant reads at SV loci → allele counts
awk 'BEGIN{OFS="\t"} !/^#/ { print $1, $2-1, $2 }' \
    $snp_dir/het_near_sv_${samp_name}.vcf \
    > $snp_dir/pos_${samp_name}.bed

$samtools view -f 1 -F 2 -b $tbam > $snp_dir/sup_${samp_name}.bam
$samtools index $snp_dir/sup_${samp_name}.bam

bcftools mpileup \
    -f $ref -a DP -A \
    -R $snp_dir/pos_${samp_name}.bed \
    $snp_dir/sup_${samp_name}.bam \
    -Ov > $snp_dir/het_on_sv_${samp_name}.vcf

rm $snp_dir/pos_${samp_name}.bed \
   $snp_dir/sup_${samp_name}.bam \
   $snp_dir/sup_${samp_name}.bam.bai

# Flat copies read by R analysis script
cp $snp_dir/het_near_sv_${samp_name}.vcf $work_dir/SNP/het_near_sv_${samp_name}.vcf
cp $snp_dir/het_on_sv_${samp_name}.vcf   $work_dir/SNP/het_on_sv_${samp_name}.vcf

# Archive raw interval VCFs
mkdir -p $snp_dir/int_snps
mv $snp_dir/int_*_SNP.vcf.gz* \
   $snp_dir/int_*_het_snp.vcf.gz* \
   $snp_dir/int_snps/ 2>/dev/null

conda deactivate
echo "Done: SNP processing [$rep / $samp_name]"
