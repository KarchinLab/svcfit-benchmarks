#!/usr/bin/env bash
#SBATCH --job-name=svclone_existing
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:15:00
#SBATCH --array=0-2249%25

# Stage immutable inputs and normalized VCF copies in a fresh result root, then
# run SVclone. This worker deliberately does not invoke Manta or SVtyper.

set -euo pipefail

script_dir=${VISOR_PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-}}
if [[ -z "$script_dir" || ! -r "$script_dir/normalize_two_sample_vcf.awk" ]]; then
    script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fi
: "${SOURCE_REPLICATES_DIR:?SOURCE_REPLICATES_DIR is required}"
: "${VISOR_REPLICATES_DIR:?VISOR_REPLICATES_DIR is required}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

source_root=${SOURCE_REPLICATES_DIR%/}
result_root=${VISOR_REPLICATES_DIR%/}

[[ "$source_root" = /* && "$result_root" = /* ]] || {
    echo "ERROR: source and result roots must be absolute paths" >&2
    exit 2
}
[[ "$source_root" != "$result_root" ]] || {
    echo "ERROR: result root must differ from the immutable source root" >&2
    exit 2
}
[[ "$SLURM_ARRAY_TASK_ID" =~ ^[0-9]+$ ]] && (( SLURM_ARRAY_TASK_ID < 2250 )) || {
    echo "ERROR: task ID must be in 0-2249: $SLURM_ARRAY_TASK_ID" >&2
    exit 2
}

rep_idx=$(( SLURM_ARRAY_TASK_ID / 75 ))
remainder=$(( SLURM_ARRAY_TASK_ID % 75 ))
exp_idx=$(( remainder / 15 ))
cond_idx=$(( remainder % 15 ))
rep=rep$(( rep_idx + 1 ))
exp_name=exp$(( exp_idx + 1 ))
purities=(10 10 10 20 20 20 40 40 40 60 60 60 80 80 80)
mixtures=(10 30 50 10 30 50 10 30 50 10 30 50 10 30 50)
cond_name=c50p${purities[$cond_idx]}m${mixtures[$cond_idx]}

source_work="$source_root/$rep/$exp_name/$cond_name"
result_work="$result_root/$rep/$exp_name/$cond_name"
source_manta="$source_work/manta/results/variants/${cond_name}.vcf"
source_svtyp="$source_work/svtyp/svt_${cond_name}.vcf"

for path in \
    "$source_manta" \
    "$source_svtyp" \
    "$source_work/short/short.out/sim.srt.bam" \
    "$source_work/facet/${cond_name}.RData" \
    "$source_work/facet/${cond_name}.bed"; do
    [[ -s "$path" ]] || { echo "ERROR: required source is missing or empty: $path" >&2; exit 3; }
done
for path in "$source_work/SNP" "$source_work/facet" "$source_work/short"; do
    [[ -d "$path" ]] || { echo "ERROR: required source directory is missing: $path" >&2; exit 3; }
done

mkdir -p "$result_work/manta/results/variants" "$result_work/svtyp"

# New runs record their complete source mapping in INPUTS.tsv and resolve
# immutable inputs directly from source_root. Result roots contain outputs,
# normalized VCF copies, and provenance only; no directory views are created.
if [[ -n "${VISOR_INPUT_MANIFEST:-}" ]]; then
    [[ -s "$VISOR_INPUT_MANIFEST" ]] || {
        echo "ERROR: VISOR_INPUT_MANIFEST is missing or empty: $VISOR_INPUT_MANIFEST" >&2
        exit 4
    }
    manifest_mapping=$(awk -F '\t' -v id="$SLURM_ARRAY_TASK_ID" '
        NR == 1 {
            for (i = 1; i <= NF; i++) column[$i] = i
            if (!("task_id" in column) || !("source_work" in column) ||
                !("result_work" in column)) exit 2
            next
        }
        $(column["task_id"]) == id {
            print $(column["source_work"]) "\t" $(column["result_work"])
            found++
        }
        END { if (found != 1) exit 1 }
    ' "$VISOR_INPUT_MANIFEST") || {
        echo "ERROR: manifest must contain exactly one row for task $SLURM_ARRAY_TASK_ID" >&2
        exit 4
    }
    IFS=$'\t' read -r manifest_source manifest_result <<< "$manifest_mapping"
    [[ "$manifest_source" == "$source_work" && "$manifest_result" == "$result_work" ]] || {
        echo "ERROR: manifest mapping does not match decoded task paths" >&2
        exit 4
    }
fi

rdata="$result_work/svclone/$cond_name/ccube_out/${cond_name}_ccube_sv_results.RData"
if [[ -s "$rdata" && "${FORCE_DOWNSTREAM_RERUN:-0}" != 1 ]]; then
    echo "Already complete: $rdata"
    exit 0
fi

normalizer="$script_dir/normalize_two_sample_vcf.awk"
[[ -r "$normalizer" ]] || { echo "ERROR: missing normalizer: $normalizer" >&2; exit 5; }
normal_sample=${rep}_${exp_name}_${cond_name}_normal
tumor_sample=${rep}_${exp_name}_${cond_name}_tumor

normalize_copy() {
    local source=$1 target=$2 tmp
    tmp="${target}.tmp.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID}"
    trap 'rm -f -- "$tmp"' RETURN
    awk -v normal_sample="$normal_sample" -v tumor_sample="$tumor_sample" \
        -f "$normalizer" "$source" > "$tmp"
    mv "$tmp" "$target"
    trap - RETURN
}

normalize_copy "$source_manta" "$result_work/manta/results/variants/${cond_name}.vcf"
normalize_copy "$source_svtyp" "$result_work/svtyp/svt_${cond_name}.vcf"

export VISOR_REPLICATES_DIR="$result_root"
export VISOR_INPUT_REPLICATES_DIR="$source_root"
bash "$script_dir/04_svclone.sh"
echo "Validated downstream result: $rdata"
