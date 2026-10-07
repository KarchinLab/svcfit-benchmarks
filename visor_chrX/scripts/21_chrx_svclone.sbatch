#!/bin/bash
#SBATCH --job-name=chrx_svclone
#SBATCH --output=log/svclone_%A_%a.out
#SBATCH --error=log/svclone_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=0-44%45
#
# SVclone across the chrX condition grid, one task per tumour BAM. The counterpart to
# 12 (calling) and 17 (SVCFit join + score); this is the competing method.
#
# WHY THIS EXISTS RATHER THAN visor_replicates/04_svclone.sh. That script is for the autosomal
# benchmark and assumes its layout (${rep}/${exp}/${cond}/manta|facet|svclone) and FACETS input.
# The chrX simulation has neither. It also sources pipeline_common.sh, which expects the
# autosomal replicate layout, and calls a make_input.R that historically assumed an unavailable
# R package. This implementation parses its arguments without
# that cross-environment dependency.
#
# FAIRNESS, since this is a head-to-head. Both methods get:
#   - the SAME call set, from the SAME file: calls*/<name>/svtyp/svt_<name>.vcf. No PASS filter,
#     because 13_chrx_score_svcfit.R applies none; filtering here would hand SVclone 32 of 42 calls
#     in a typical condition and lose it detections it never had a chance at.
#   - the SAME copy-number evidence: the DNAcopy segments in seg10k*/, since neither method has
#     FACETS and FACETS could not work on a male X anyway.
# SVclone additionally gets TRUE PURITY, which SVCFit never receives (calc_svcf takes no purity
# argument). That asymmetry favours SVclone and matches the published autosomal benchmark. It must
# be stated wherever these results are reported.
#
# Needs COV and, for a replicate, REP in the environment; both come through --export.
# Optional: SVC_FLOOR (default 1) for the copy-number floor sensitivity run.

set -euo pipefail

_visor_walk_up() {
    local d="${SLURM_SUBMIT_DIR:-$PWD}"
    while [[ "$d" != "/" && -n "$d" ]]; do
        [[ -r "$d/config.local.sh" ]] && { printf '%s' "$d/visor_chrX/scripts"; return 0; }
        d="$(dirname "$d")"
    done
    return 1
}
_resolve_scripts_dir() {
    local d
    for d in "${CHRX_SCRIPTS:-}" "${SLURM_SUBMIT_DIR:-}/scripts" "${SLURM_SUBMIT_DIR:-}" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" "$(_visor_walk_up)"; do
        [[ -n "$d" && -r "$d/chrx_common.sh" ]] && { printf '%s' "$d"; return 0; }
    done
    return 1
}
SCRIPTS="$(_resolve_scripts_dir)" || { echo "ERROR: cannot locate scripts dir" >&2; exit 1; }
source "${SCRIPTS}/chrx_common.sh"          # COV / REP -> cov, rep_tag, BASE_DIR, cond_tsv

V="$BASE_DIR"
tid="${SLURM_ARRAY_TASK_ID:-${1:-}}"
[[ -n "$tid" ]] || { echo "ERROR: no task id (pass one as \$1 outside SLURM)" >&2; exit 1; }

conds="$cond_tsv"
[[ -s "$conds" ]] || { echo "ERROR: missing $conds" >&2; exit 1; }

# The conditions table is the single source for exp, condition, TRUE purity and the tumour BAM, so
# the four cannot drift apart the way they would if purity were re-parsed out of the condition name.
read -r _tid exp cond ppur _mix _c2 _c3 _norm _seed tbam < <(awk -F'\t' -v t="$tid" 'NR>1 && $1==t' "$conds")
[[ -n "${exp:-}" && -n "${cond:-}" ]] || { echo "ERROR: no condition for task $tid in $conds" >&2; exit 1; }
name="${exp}_${cond}"

vcf="${V}/calls${rep_tag}/${name}/svtyp/svt_${name}.vcf"
seg="${V}/seg10k${rep_tag}/chrx_seg_${name}.csv"
out="${V}/svclone${rep_tag}/${name}"
cfg="${SVCLONE_CFG:-${V}/resources/svclone_config_chrx.ini}"
floor="${SVC_FLOOR:-1}"

for f in "$vcf" "$seg" "$tbam" "$cfg"; do
    [[ -s "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

# Guard the rm -rf: only ever inside our own tree, never a bare or parent path.
if [[ -z "$out" || "$out" != "${V}/svclone"* ]]; then
    echo "ERROR: refusing to clear '$out' -- not under ${V}/svclone" >&2; exit 1
fi
rm -rf "$out"; mkdir -p "$out"; cd "$out"

echo "========================================================"
echo "  task    : ${tid}   replicate: ${REP:-0}"
echo "  cond    : ${name}   purity: ${ppur}%   floor: ${floor}"
echo "  vcf     : ${vcf}"
echo "  seg     : ${seg}"
echo "  bam     : ${tbam}"
echo "========================================================"

# SVclone lives in its own project environment and is activated by absolute prefix.
# shellcheck disable=SC1091
source "$CONDA_SH"
conda activate "$ENV_SVCLONE"
command -v svclone >/dev/null || { echo "ERROR: svclone not on PATH after activate" >&2; exit 1; }
echo "  svclone : $(svclone --version 2>&1 | head -1)"

ppur_frac=$(awk -v p="$ppur" 'BEGIN{printf "%.3f", p/100}')

echo "--- building inputs (DNAcopy -> SVclone) ---"
Rscript "${SCRIPTS}/20_chrx_svclone_input.R" \
    --vcf "$vcf" --seg "$seg" --sample "$name" \
    --purity "$ppur_frac" --out "$out" --floor "$floor"

echo "--- annotate ---"
svclone annotate -i "svc_${name}_simple.txt" -b "$tbam" -s "$name" \
    --sv_format simple -cfg "$cfg"

echo "--- count ---"
svclone count -i "${name}/${name}_svin.txt" -b "$tbam" -s "$name" -cfg "$cfg"

echo "--- filter ---"
svclone filter -s "$name" -i "${name}/${name}_svinfo.txt" \
    -p "spp_${name}.txt" -c cnv.txt -cfg "$cfg"

echo "--- cluster ---"
svclone cluster -s "$name" -cfg "$cfg"

rdata="${name}/ccube_out/${name}_ccube_sv_results.RData"
[[ -f "$rdata" ]] || { echo "ERROR: $rdata not produced" >&2; exit 1; }
echo "RESULT: SVclone ${name} (replicate ${REP:-0}) -> ${out}/${rdata}"
