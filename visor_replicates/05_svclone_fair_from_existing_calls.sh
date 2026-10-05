#!/usr/bin/env bash
#SBATCH --job-name=svclone_fair
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:15:00
#SBATCH --array=0-2249%25

# FACETS-purity SVclone arm. Existing calls are copied and relabelled; no caller
# is invoked. Missing FACETS purity is a recorded, successful ledger outcome.
set -euo pipefail

: "${SOURCE_REPLICATES_DIR:?}" "${FAIR_RESULT_ROOT:?}" "${VISOR_PIPELINE_DIR:?}"
: "${FACETS_PURITY_TABLE:?}" "${VISOR_SVCLONE_CONDA_SH:?}" "${VISOR_SVCLONE_ENV:?}"
: "${SLURM_ARRAY_TASK_ID:?}"
source_root=${SOURCE_REPLICATES_DIR%/}; result_root=${FAIR_RESULT_ROOT%/}; pipe_dir=${VISOR_PIPELINE_DIR%/}
[[ "$source_root" = /* && "$result_root" = /* && "$source_root" != "$result_root" ]] || exit 2
[[ "$SLURM_ARRAY_TASK_ID" =~ ^[0-9]+$ ]] && ((SLURM_ARRAY_TASK_ID < 2250)) || exit 2

rep_idx=$((SLURM_ARRAY_TASK_ID/75)); rem=$((SLURM_ARRAY_TASK_ID%75)); exp_idx=$((rem/15)); ci=$((rem%15))
pur=(10 10 10 20 20 20 40 40 40 60 60 60 80 80 80)
mix=(10 30 50 10 30 50 10 30 50 10 30 50 10 30 50)
rep=rep$((rep_idx+1)); exp=exp$((exp_idx+1)); cond=c50p${pur[$ci]}m${mix[$ci]}
src="$source_root/$rep/$exp/$cond"; dst="$result_root/$rep/$exp/$cond"; svc="$dst/svclone_fair"
vcf_src="$src/manta/results/variants/$cond.vcf"; vcf="$dst/manta/results/variants/$cond.vcf"
for f in "$vcf_src" "$src/short/short.out/sim.srt.bam" "$src/facet/$cond.bed"; do
  [[ -s "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 3; }
done

rdata="$svc/$cond/ccube_out/${cond}_ccube_sv_results.RData"; status="$svc/facets_status_${cond}.txt"
if [[ "${FORCE_FAIR_RERUN:-0}" != 1 ]]; then
  [[ -s "$rdata" ]] && { echo "Already complete: $rdata"; exit 0; }
  [[ -s "$status" ]] && grep -q $'failed_no_purity' "$status" && { echo "Already recorded: $status"; exit 0; }
fi
[[ -z "$svc" || "$svc" != "$dst"/* ]] && { echo "ERROR: unsafe output: $svc" >&2; exit 4; }
rm -rf -- "$svc"; mkdir -p "$svc" "$(dirname "$vcf")"

normal="$pipe_dir/normalize_two_sample_vcf.awk"; fair_parser="$pipe_dir/helpers/make_input_fair.R"
tmp="${vcf}.tmp.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID}"; trap 'rm -f -- "$tmp"' EXIT
awk -v normal_sample="${rep}_${exp}_${cond}_normal" -v tumor_sample="${rep}_${exp}_${cond}_tumor" \
  -f "$normal" "$vcf_src" > "$tmp"
mv "$tmp" "$vcf"; trap - EXIT

export PATH="$VISOR_SVCLONE_ENV/bin:$PATH"
rscript="$VISOR_SVCLONE_ENV/bin/Rscript"; svclone="$VISOR_SVCLONE_ENV/bin/svclone"
[[ -x "$rscript" && -x "$svclone" ]] || { echo "ERROR: incomplete SVclone environment: $VISOR_SVCLONE_ENV" >&2; exit 4; }
set +e
"$rscript" "$fair_parser" -f "$vcf" -p "0.${pur[$ci]}" -o "$svc" -c "$src/facet" \
  -s "$cond" -r "$rep" -e "$exp" -t "$FACETS_PURITY_TABLE"
rc=$?
set -e
if ((rc == 3)); then echo "Done: expected FACETS-purity failure [$rep/$exp/$cond]"; exit 0; fi
((rc == 0)) || { echo "ERROR: fair parser rc=$rc" >&2; exit "$rc"; }

cd "$svc"; cfg="$pipe_dir/helpers/svclone_config.ini"; bam="$src/short/short.out/sim.srt.bam"
"$svclone" annotate -i "svc_${cond}_simple.txt" -b "$bam" -s "$cond" --sv_format simple -cfg "$cfg"
"$svclone" count -i "$cond/${cond}_svin.txt" -b "$bam" -s "$cond" -cfg "$cfg"
"$svclone" filter -s "$cond" -i "$cond/${cond}_svinfo.txt" -p "spp_${cond}.txt" -c cnv.txt -cfg "$cfg"
"$svclone" cluster -s "$cond" -cfg "$cfg"
"$rscript" -e 'load(commandArgs(TRUE)[1]); stopifnot(exists("doubleBreakPtsRes"))' "$rdata"
echo "Done: SVclone FAIR [$rep/$exp/$cond]"
