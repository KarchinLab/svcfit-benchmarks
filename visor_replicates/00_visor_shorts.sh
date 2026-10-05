#!/bin/bash
#SBATCH --job-name=visor_rep
#SBATCH --output=log/00_visor_%A_%a.out
#SBATCH --error=log/00_visor_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=15G
#SBATCH --time=8:00:00
#SBATCH --array=0-2249%10
#
# Generate VISOR SHORtS tumor BAMs with replicate-specific seeds.
# 2250 tasks = 30 replicates x 5 experiments x 15 conditions
#             (p10/p20/p40/p60/p80 x m10/m30/m50).
# HACk genomes and normal BAMs are reused (only read noise varies).
# NOTE: All 2250 BAMs are already generated. Do not re-run unless
#       starting completely from scratch.

set -euo

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
set +o pipefail   # VISOR prints non-fatal warnings that would trip pipefail
log_repro_header "00_visor_shorts.sh"

source "$CONDA_SH"
conda activate "$ENV_VISOR"
echo "  VISOR   : $(VISOR --version 2>&1 | head -1)"
echo "  samtools: $($samtools --version | head -1)"

threads=$SLURM_CPUS_PER_TASK

mkdir -p $short_dir
rm -rf $short_dir/short.bed $short_dir/short.out

# Build BED file for VISOR SHORtS (max coordinates across all sources + reference)
cut -f1,2 \
    $hack_dir/c2/*.fai $hack_dir/c3/*.fai $snp_hack/*.fai $ref.fai \
    | sort \
    | awk '$2 > maxvals[$1] {lines[$1]=$0; maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' \
    | awk -v p="${pur}" 'OFS=FS="\t"''{print $1, "1", $2, '100', p}' > $short_dir/short.bed

echo "Running VISOR SHORtS with seed=${visor_seed}..."

"$ENV_VISOR/bin/python" $seeded_py $visor_seed \
    -g $ref \
    -s $hack_dir/c2 $hack_dir/c3 $snp_hack \
    -b $short_dir/short.bed \
    -o $short_dir/short.out \
    --mutation 0 \
    --extindels 0 \
    --threads $threads \
    --coverage $cov \
    --clonefraction ${frac[@]} \
    --error 0 \
    --indels 0

if [[ ! -f $tbam ]]; then
    echo "ERROR: $tbam not produced" >&2; exit 1
fi

echo "Done: $tbam ($(du -sh $tbam 2>/dev/null | cut -f1))"
conda deactivate
