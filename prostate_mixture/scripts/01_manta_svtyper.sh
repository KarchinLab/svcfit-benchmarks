#!/bin/bash
#SBATCH --job-name=prst_manta
#SBATCH --output=log/01_manta_%A_%a.out
#SBATCH --error=log/01_manta_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --array=0-329%8

set -euo pipefail

# Locate pipeline_common.sh beside this script (same lookup as 04a/04b).
_pc=""
for _c in "${SLURM_SUBMIT_DIR:-}" "${PWD:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
    [[ -n "$_c" && -r "$_c/pipeline_common.sh" ]] && { _pc="$_c/pipeline_common.sh"; break; }
done
[[ -n "$_pc" ]] || { echo "ERROR: cannot locate pipeline_common.sh" >&2; exit 1; }
# shellcheck disable=SC1090
source "$_pc"
mkdir -p $manta_dir $svtyp_dir $work_dir/prostate_vcfs

if [[ -f "$work_dir/prostate_vcfs/svt_${samp_name}.vcf" ]]; then
    echo "[$rep / $samp_name] Output already exists, skipping: svt_${samp_name}.vcf"
    exit 0
fi

log_repro_header "01_manta_svtyper.sh"

echo "[$rep / $samp_name] Manta + SVtyper"

# ---- Manta ----
source "$PIPELINE_CONDA_SH"
set +u  # Manta activation references unset JAVA_HOME/JAVA_LD_LIBRARY_PATH variables.
conda activate "$ENV_MANTA"
set -u
[[ -x "$MANTA_CONFIG" && -x "$MANTA_CONVERT" ]] || {
    echo "ERROR: project Manta executables are missing" >&2; exit 1;
}
echo "  manta   : $($ENV_MANTA/bin/python "$MANTA_CONFIG" --version 2>&1 | head -1)"

# Remove any incomplete manta run directory; configManta.py fails if it already exists.
# Safe here because the skip guard above already confirmed the final output is absent.
rm -rf "$manta_dir"
mkdir -p "$manta_dir"

"$ENV_MANTA/bin/python" "$MANTA_CONFIG" \
    --normalBam $nbam \
    --tumorBam  $tbam \
    --referenceFasta $ref \
    --runDir $manta_dir

"$ENV_MANTA/bin/python" $manta_dir/runWorkflow.py -j 6

cd $manta_dir/results/variants/
"$ENV_MANTA/bin/python" "$MANTA_CONVERT" \
    $samtools $ref somaticSV.vcf.gz > ${samp_name}.vcf

conda deactivate

# ---- SVtyper ----
source "$PIPELINE_CONDA_SH"
conda activate "$ENV_SVTYPER"
echo "  svtyper : $(svtyper --version 2>&1 | head -1)"

awk 'BEGIN{FS=OFS="\t"}
     /^#/ { print; next }
     { info=$8
       gsub(/CIPOS=[^;]+;?/,"",info)
       gsub(/CIEND=[^;]+;?/,"",info)
       $8="CIPOS=-100,100;CIEND=-100,100;"info
       print }' \
    $manta_dir/results/variants/${samp_name}.vcf \
    > $svtyp_dir/${samp_name}_w_CI.vcf

svtyper \
    -i $svtyp_dir/${samp_name}_w_CI.vcf \
    -B $tbam \
    -l $svtyp_dir/${samp_name}.bam.json \
    > $svtyp_dir/svt_${samp_name}.vcf

# Flat copy read by R analysis script
cp $svtyp_dir/svt_${samp_name}.vcf $work_dir/prostate_vcfs/svt_${samp_name}.vcf

conda deactivate
echo "Done: Manta + SVtyper [$rep / $samp_name]"
