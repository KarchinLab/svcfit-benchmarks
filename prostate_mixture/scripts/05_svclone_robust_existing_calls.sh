#!/usr/bin/env bash
#SBATCH --job-name=prst_svclone_fix
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:15:00
#SBATCH --array=0-329%20
set -euo pipefail
: "${PROSTATE_SOURCE_ROOT:?}" "${PROSTATE_MIX_BAM_ROOT:?}" "${PROSTATE_RESULT_ROOT:?}"
: "${PROSTATE_PIPELINE_DIR:?}" "${VISOR_SVCLONE_CONDA_SH:?}" "${VISOR_SVCLONE_ENV:?}" "${SLURM_ARRAY_TASK_ID:?}"
conditions=(3m19 3m28 3m37 3m46 3m55 3m64 3m73 3m82 3m91 4m 5m)
rep=rep$((SLURM_ARRAY_TASK_ID/11+1)); cond=${conditions[$((SLURM_ARRAY_TASK_ID%11))]}
src="${PROSTATE_SOURCE_ROOT%/}/$rep"; dst="${PROSTATE_RESULT_ROOT%/}/$rep"; svc="$dst/svclone/$cond"
vcf_src="$src/manta/$cond/results/variants/$cond.vcf"; facet="$src/facet/$cond"; bam="${PROSTATE_MIX_BAM_ROOT%/}/$rep/$cond.bam"
for f in "$vcf_src" "$facet/$cond.RData" "$facet/$cond.bed" "$bam";do [[ -s "$f" ]]||{ echo "ERROR: missing $f" >&2;exit 3;};done
rdata="$svc/$cond/ccube_out/${cond}_ccube_sv_results.RData"
[[ -s "$rdata" && "${FORCE_PROSTATE_RERUN:-0}" != 1 ]]&&{ echo "Already complete: $rdata";exit 0;}
[[ -z "$svc"||"$svc" != "$dst"/* ]]&&exit 4
rm -rf -- "$svc";mkdir -p "$svc" "$dst/manta/$cond/results/variants"
normal="${PROSTATE_PIPELINE_DIR%/}/../../visor_replicates/normalize_two_sample_vcf.awk"
vcf="$dst/manta/$cond/results/variants/$cond.vcf";tmp="${vcf}.tmp.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID}";trap 'rm -f -- "$tmp"' EXIT
awk -v normal_sample="${rep}_${cond}_normal" -v tumor_sample="${rep}_${cond}_tumor" -f "$normal" "$vcf_src">"$tmp";mv "$tmp" "$vcf";trap - EXIT
export PATH="$VISOR_SVCLONE_ENV/bin:$PATH"
rscript="$VISOR_SVCLONE_ENV/bin/Rscript";svclone="$VISOR_SVCLONE_ENV/bin/svclone"
[[ -x "$rscript"&&-x "$svclone" ]]||{ echo "ERROR: incomplete SVclone environment: $VISOR_SVCLONE_ENV" >&2;exit 4;}
parser="${PROSTATE_PIPELINE_DIR%/}/helper/make_input_robust.R";cfg="${PROSTATE_PIPELINE_DIR%/}/helper/svclone_config.ini"
"$rscript" "$parser" -f "$vcf" -o "$svc" -c "$facet" -s "$cond"
cd "$svc"
"$svclone" annotate -i "svc_${cond}_simple.txt" -b "$bam" -s "$cond" --sv_format simple -cfg "$cfg"
"$svclone" count -i "$cond/${cond}_svin.txt" -b "$bam" -s "$cond" -cfg "$cfg"
"$svclone" filter -s "$cond" -i "$cond/${cond}_svinfo.txt" -p "spp_${cond}.txt" -c cnv.txt -cfg "$cfg"
"$svclone" cluster -s "$cond" -cfg "$cfg"
"$rscript" -e 'load(commandArgs(TRUE)[1]);stopifnot(exists("doubleBreakPtsRes"))' "$rdata"
echo "Done: prostate SVclone robust [$rep/$cond]"
