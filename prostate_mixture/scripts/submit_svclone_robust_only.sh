#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")"&&pwd);src="";bam="";out="";conda_sh="";env="";array="0,9,10%3";submit=0
while(($#));do case "$1" in --source-root)src=${2:?};shift 2;;--bam-root)bam=${2:?};shift 2;;--result-root)out=${2:?};shift 2;;--svclone-conda-sh)conda_sh=${2:?};shift 2;;--svclone-env)env=${2:?};shift 2;;--array)array=${2:?};shift 2;;--submit)submit=1;shift;;*)exit 2;;esac;done
for p in "$src" "$bam" "$out" "$conda_sh" "$env";do [[ "$p" = /* ]]||exit 2;done
[[ -d "$src"&&-d "$bam"&&-r "$conda_sh"&&-d "$env"&&"$src" != "$out" ]]||exit 2
missing=0;conds=(3m19 3m28 3m37 3m46 3m55 3m64 3m73 3m82 3m91 4m 5m)
for r in $(seq 1 30);do for c in "${conds[@]}";do for f in "$src/rep$r/manta/$c/results/variants/$c.vcf" "$src/rep$r/facet/$c/$c.RData" "$src/rep$r/facet/$c/$c.bed" "$bam/rep$r/$c.bam";do [[ -s "$f" ]]||missing=$((missing+1));done;done;done
((missing==0))||{ echo "ERROR: missing inputs=$missing" >&2;exit 3;}
echo "Preflight passed: 330 conditions; calling disabled; array=$array; result=$out"
((submit))||{ echo "DRY RUN ONLY";exit 0;};mkdir -p "$out/log"
job=$(sbatch --parsable --array="$array" --output="$out/log/svclone_%A_%a.out" --error="$out/log/svclone_%A_%a.err" --export="ALL,PROSTATE_SOURCE_ROOT=$src,PROSTATE_MIX_BAM_ROOT=$bam,PROSTATE_RESULT_ROOT=$out,PROSTATE_PIPELINE_DIR=$here,VISOR_SVCLONE_CONDA_SH=$conda_sh,VISOR_SVCLONE_ENV=$env" "$here/05_svclone_robust_existing_calls.sh")
printf 'job=%s\ncommit=%s\narray=%s\ncalling=disabled\n' "$job" "$(git -C "$here" rev-parse HEAD)" "$array">"$out/SUBMISSION-${job}.txt";echo "Submitted $job"
