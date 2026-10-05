#!/bin/bash
#SBATCH --job-name=manta_vis
#SBATCH --output=log/01_manta_%A_%a.out
#SBATCH --error=log/01_manta_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=12G
#SBATCH --time=8:00:00
#SBATCH --array=0-2249%10
#
# Manta SV calling + SVtyper genotyping on VISOR replicate BAMs.

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
log_repro_header "01_manta_svtyper.sh"
normalizer="$(dirname "$_pc")/normalize_two_sample_vcf.awk"
[[ -r "$normalizer" ]] || {
    echo "ERROR: VCF sample normalizer not found: $normalizer" >&2
    exit 1
}


threads=$SLURM_CPUS_PER_TASK
mkdir -p $manta_dir $svtyp_dir

# ---- Manta ----
[[ -r "$CONDA_SH" ]] || { echo "ERROR: CONDA_SH is not readable: $CONDA_SH" >&2; exit 1; }
# shellcheck disable=SC1090
set +u
source "$CONDA_SH"
conda activate "$ENV_MANTA"
set -u
manta_bin="$MANTA_CONFIG"
manta_python="$ENV_MANTA/bin/python"
[[ -x "$manta_bin" && -x "$MANTA_CONVERT" && -x "$manta_python" ]] || {
    echo "ERROR: configured Manta tools are missing" >&2; exit 1;
}
echo "  manta   : $("$manta_python" "$manta_bin" --version 2>&1 | head -1)"

"$manta_python" "$manta_bin" \
    --normalBam $nbam \
    --tumorBam $tbam \
    --referenceFasta $ref \
    --runDir $manta_dir

"$manta_python" "$manta_dir/runWorkflow.py" -j "$threads"

converted_vcf="$manta_dir/results/variants/${cond_name}.vcf"
converted_tmp="${converted_vcf}.tmp.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"
normal_sample="${rep}_${exp_name}_${cond_name}_normal"
tumor_sample="${rep}_${exp_name}_${cond_name}_tumor"
trap 'rm -f -- "$converted_tmp"' EXIT

"$manta_python" "$MANTA_CONVERT" \
    "$samtools" "$ref" \
    "$manta_dir/results/variants/somaticSV.vcf.gz" \
    | awk -v normal_sample="$normal_sample" -v tumor_sample="$tumor_sample" \
        -f "$normalizer" \
    > "$converted_tmp"
mv "$converted_tmp" "$converted_vcf"
trap - EXIT

conda deactivate

# ---- SVtyper ----
conda activate "$ENV_SVTYPER"
echo "  svtyper : $(svtyper --version 2>&1 | head -1)"

awk 'BEGIN{FS=OFS="\t"}
     /^#/ { print; next }
     { info=$8
       gsub(/CIPOS=[^;]+;?/,"",info)
       gsub(/CIEND=[^;]+;?/,"",info)
       $8="CIPOS=-100,100;CIEND=-100,100;"info
       print }' \
    "$converted_vcf" \
    > "$svtyp_dir/${cond_name}_w_CI.vcf"

ci_vcf="$svtyp_dir/${cond_name}_w_CI.vcf"
final_vcf="$svtyp_dir/svt_${cond_name}.vcf"
compat_vcf="${ci_vcf}.svtyper-input.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"
raw_vcf="${final_vcf}.raw.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"
final_tmp="${final_vcf}.tmp.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"

# SVtyper identifies its target sample from the BAM read-group SM value. The
# immutable VISOR BAMs use SM:bulk, so a pre-normalized VCF would otherwise
# make SVtyper append a third sample named "bulk". Use the BAM name only in a
# temporary compatibility VCF, then restore the deterministic public names.
bam_sample="$("$samtools" view -H "$tbam" | awk -F '\t' '
    $1 == "@RG" {
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^SM:/) samples[substr($i, 4)] = 1
        }
    }
    END {
        for (sample in samples) { count++; only = sample }
        if (count != 1) exit 1
        print only
    }')" || {
    echo "ERROR: tumor BAM must contain exactly one read-group SM value: $tbam" >&2
    exit 1
}

trap 'rm -f -- "$compat_vcf" "$raw_vcf" "$final_tmp"' EXIT
awk -v normal_sample="$normal_sample" -v tumor_sample="$bam_sample" \
    -v emit_metadata=0 -f "$normalizer" "$ci_vcf" > "$compat_vcf"

svtyper \
    -i "$compat_vcf" \
    -B "$tbam" \
    -l "$svtyp_dir/${cond_name}.bam.json" \
    > "$raw_vcf"

awk -v normal_sample="$normal_sample" -v tumor_sample="$tumor_sample" \
    -v emit_metadata=0 -f "$normalizer" "$raw_vcf" > "$final_tmp"
mv "$final_tmp" "$final_vcf"
rm -f -- "$compat_vcf" "$raw_vcf"
trap - EXIT

conda deactivate
echo "Done: Manta + SVtyper [$rep / $exp_name / $cond_name]"
