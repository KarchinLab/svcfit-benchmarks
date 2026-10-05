#!/bin/bash
#SBATCH --job-name=make_mix
#SBATCH --output=log/00_make_mix_%A_%a.out
#SBATCH --error=log/00_make_mix_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --array=58-59,63,65-84,86-87,94,96,98-99,101-107
####0-329
set -euo pipefail
# Creates 30 independent replicate in-silico mixture BAMs from the original
# unmixed BAMs, each using a different random seed for samtools subsampling.
#
# Seed scheme (systematic, no overlap with original published data):
#   3-clone  : single seed per replicate = rep_idx * 111 + 11
#              (rep1=11, rep2=122, rep3=233, rep4=344, rep5=455)
#   4-clone  : three seeds = (rep_idx+1)*100 + {1,2,3}
#   5-clone  : four seeds  = (rep_idx+1)*100 + {11,12,13,14}
#
# samtools -s format: SEED.FRACTION
#   e.g. -s 11.3  → seed=11, keep 30% of reads

# Locate pipeline_common.sh beside this script (same lookup as 04a/04b).
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"

# $work_dir is defined once, in pipeline_common.sh.
out_dir=$work_dir
mkdir -p $out_dir

threads=8
mem=1   # per-thread sort buffer; 9 threads × 1G = 9G total, within --mem=12G

if [[ -f "$tbam" && -f "${tbam}.bai" ]]; then
    echo "[$rep / $samp_name] Output already exists, skipping: $tbam"
    exit 0
fi

log_repro_header "00_make_mixtures.sh"
echo "  samtools: $($samtools --version | head -1)"
echo "  bM_bam size : $(du -sh $bM_bam 2>/dev/null | cut -f1)"
echo "  gM_bam size : $(du -sh $gM_bam 2>/dev/null | cut -f1)"
echo ""

source "$PIPELINE_CONDA_SH"
conda activate "$ENV_VISOR"

echo "[$rep / $samp_name] Creating mixture BAM"

# ============================================================
# 3-clone mixtures (conditions 0-8)
# ============================================================
if [ $cond_idx -le 8 ]; then
    seed=$(( rep_idx * 111 + 11 ))
    bm_digit=$(echo $samp_name | sed 's/3m\(.\)./\1/')  # 1st digit = bM%/10
    gm_digit=$(echo $samp_name | sed 's/3m.\(.\)/\1/')  # 2nd digit = gM%/10

    echo "  seed=$seed  bM=${bm_digit}0%  gM=${gm_digit}0%"

    $samtools view -@ $threads -s ${seed}.${bm_digit} -b $bM_bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_bM_sub.bam
    $samtools view -@ $threads -s ${seed}.${gm_digit} -b $gM_bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_gM_sub.bam

    $samtools merge -n -@ $threads ${out_dir}/${samp_name}_merged.bam \
        ${out_dir}/${samp_name}_bM_sub.bam ${out_dir}/${samp_name}_gM_sub.bam
    $samtools sort  -@ $threads -m ${mem}G ${out_dir}/${samp_name}_merged.bam -o ${out_dir}/${samp_name}_sorted.bam
    rm ${out_dir}/${samp_name}_bM_sub.bam ${out_dir}/${samp_name}_gM_sub.bam ${out_dir}/${samp_name}_merged.bam
    # Unify SM tags so bcftools mpileup outputs one sample column
    $samtools view -H ${out_dir}/${samp_name}_sorted.bam \
        | awk '/^@RG/{gsub(/SM:[^\t]+/, "SM:tumor")} {print}' \
        | $samtools reheader - ${out_dir}/${samp_name}_sorted.bam > $tbam
    rm ${out_dir}/${samp_name}_sorted.bam
    $samtools index -@ $threads $tbam

# ============================================================
# 4-clone mixture (condition 9): 4m
# bM odd-chrom 20%, gM odd-chrom 40%, bM even-chrom 60%
# ============================================================
elif [ $cond_idx -eq 9 ]; then
    base=$(( (rep_idx + 1) * 100 ))
    sa=$(( base + 1 )); sb=$(( base + 2 )); sc=$(( base + 3 ))
    echo "  4-clone seeds: $sa $sb $sc"

    $samtools view -@ $threads -s ${sa}.20 -b ${split}/fil_bM_odd.bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_p1.bam
    $samtools view -@ $threads -s ${sb}.40 -b ${split}/fil_gM_odd.bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_p2.bam
    $samtools view -@ $threads -s ${sc}.60 -b ${split}/fil_bM_even.bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_p3.bam

    $samtools merge -n -@ $threads ${out_dir}/${samp_name}_merged.bam \
        ${out_dir}/${samp_name}_p1.bam ${out_dir}/${samp_name}_p2.bam ${out_dir}/${samp_name}_p3.bam
    $samtools sort  -@ $threads -m ${mem}G ${out_dir}/${samp_name}_merged.bam -o ${out_dir}/${samp_name}_sorted.bam
    rm ${out_dir}/${samp_name}_p1.bam ${out_dir}/${samp_name}_p2.bam ${out_dir}/${samp_name}_p3.bam ${out_dir}/${samp_name}_merged.bam
    $samtools view -H ${out_dir}/${samp_name}_sorted.bam \
        | awk '/^@RG/{gsub(/SM:[^\t]+/, "SM:tumor")} {print}' \
        | $samtools reheader - ${out_dir}/${samp_name}_sorted.bam > $tbam
    rm ${out_dir}/${samp_name}_sorted.bam
    $samtools index -@ $threads $tbam

# ============================================================
# 5-clone mixture (condition 10): 5m
# bM odd 80%, bM even 60%, gM odd 20%, gM even 40%
# ============================================================
elif [ $cond_idx -eq 10 ]; then
    base=$(( (rep_idx + 1) * 100 ))
    sa=$(( base + 11 )); sb=$(( base + 12 ))
    sc=$(( base + 13 )); sd=$(( base + 14 ))
    echo "  5-clone seeds: $sa $sb $sc $sd"

    $samtools view -@ $threads -s ${sa}.80 -b ${split}/fil_bM_odd.bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_p1.bam
    $samtools view -@ $threads -s ${sb}.60 -b ${split}/fil_bM_even.bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_p2.bam
    $samtools view -@ $threads -s ${sc}.20 -b ${split}/fil_gM_odd_downs.bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_p3.bam
    $samtools view -@ $threads -s ${sd}.40 -b ${split}/fil_gM_even_downs.bam \
        | $samtools sort -@ $threads -m ${mem}G -o ${out_dir}/${samp_name}_p4.bam

    $samtools merge -n -@ $threads ${out_dir}/${samp_name}_merged.bam \
        ${out_dir}/${samp_name}_p1.bam ${out_dir}/${samp_name}_p2.bam ${out_dir}/${samp_name}_p3.bam ${out_dir}/${samp_name}_p4.bam
    $samtools sort  -@ $threads -m ${mem}G ${out_dir}/${samp_name}_merged.bam -o ${out_dir}/${samp_name}_sorted.bam
    rm ${out_dir}/${samp_name}_p1.bam ${out_dir}/${samp_name}_p2.bam ${out_dir}/${samp_name}_p3.bam \
       ${out_dir}/${samp_name}_p4.bam ${out_dir}/${samp_name}_merged.bam
    $samtools view -H ${out_dir}/${samp_name}_sorted.bam \
        | awk '/^@RG/{gsub(/SM:[^\t]+/, "SM:tumor")} {print}' \
        | $samtools reheader - ${out_dir}/${samp_name}_sorted.bam > $tbam
    rm ${out_dir}/${samp_name}_sorted.bam
    $samtools index -@ $threads $tbam
fi

echo "Done: $tbam  ($(du -sh $tbam 2>/dev/null | cut -f1))"
conda deactivate
