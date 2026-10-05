#!/bin/bash
#SBATCH --job-name=fil_src
#SBATCH --output=log/00a_filter_%A_%a.out
#SBATCH --error=log/00a_filter_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --array=0-6
#
# Pre-filter all source BAMs to uniform read length (L=100) required by SVclone.
# SVclone's bamtools.estimateTagSize() raises ValueError on mixed read lengths.
# Run ONCE before 00_make_mixtures.sh; output saved in source_bams/ (originals untouched).
#
# Array index → BAM:
#   0  bM_recal_sorted.bam        (~130G)
#   1  gM_recal_sorted.bam        (~130G)
#   2  bM_odd_chroms.bam
#   3  gM_odd_chroms.bam
#   4  bM_even_chroms.bam
#   5  gM_odd_chroms_downs.bam
#   6  gM_even_chroms_downs.bam

set -euo pipefail

# Paths come from pipeline_common.sh, like every other stage. This runs before the
# array decoding is meaningful, so SLURM_ARRAY_TASK_ID is defaulted for sourcing.
: "${SLURM_ARRAY_TASK_ID:=0}"
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"

mkdir -p "$src_bams"
split="$split_chrom"   # the UNfiltered split-chrom BAMs are this script's input

source "$CONDA_SH"
conda activate "$ENV_VISOR"   # has samtools

declare -a SRC_LIST=(
    "${og_bam}/bM_recal_sorted.bam"
    "${og_bam}/gM_recal_sorted.bam"
    "${split}/bM_odd_chroms.bam"
    "${split}/gM_odd_chroms.bam"
    "${split}/bM_even_chroms.bam"
    "${split}/gM_odd_chroms_downs.bam"
    "${split}/gM_even_chroms_downs.bam"
)

declare -a DST_LIST=(
    "${src_bams}/fil_bM.bam"
    "${src_bams}/fil_gM.bam"
    "${src_bams}/fil_bM_odd.bam"
    "${src_bams}/fil_gM_odd.bam"
    "${src_bams}/fil_bM_even.bam"
    "${src_bams}/fil_gM_odd_downs.bam"
    "${src_bams}/fil_gM_even_downs.bam"
)

src=${SRC_LIST[$SLURM_ARRAY_TASK_ID]}
dst=${DST_LIST[$SLURM_ARRAY_TASK_ID]}

echo ""
echo "========================================================"
echo "  SCRIPT  : 00a_filter_source_bams.sh"
echo "  DATE    : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  HOST    : $(hostname -s)"
echo "  SLURM_JOB_ID        : ${SLURM_JOB_ID:-N/A}"
echo "  SLURM_ARRAY_JOB_ID  : ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "  SLURM_ARRAY_TASK_ID : ${SLURM_ARRAY_TASK_ID}"
echo "  samtools: $($samtools --version | head -1)"
echo "  filter  : length == 100 bp"
echo "  src     : ${src}  ($(du -sh $src 2>/dev/null | cut -f1))"
echo "  dst     : ${dst}"
echo "========================================================"
echo ""
echo "Filtering: $src → $dst"

# Streaming pipe: decompress → awk length filter → re-compress
# -@ 1 gives 2 threads per samtools process; awk runs concurrently as 3rd process.
$samtools view -h -@ 1 "$src" \
    | awk -v L=100 'BEGIN{OFS="\t"} /^@/{print; next} length($10)==L{print}' \
    | $samtools view -b -@ 1 -o "$dst"

$samtools index "$dst"

echo "Done: $dst"
conda deactivate
