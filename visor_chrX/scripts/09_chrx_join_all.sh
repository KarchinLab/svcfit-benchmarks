#!/bin/bash
# Run the segment-to-SV join (script 08) across the whole condition grid.
#
# For the simulation, one "sample" is one (experiment, condition) pair, and the
# SV list comes from truth/chrx_sv_ground_truth_c<cov>.tsv rather than a caller --
# Manta and SVtyper are not installed, and calling is out of scope. That means
# this exercises the JOIN, not SV detection: breakpoints are exact, so any
# cn_bar error here is the join's or the segmentation's, never a caller's.
#
# Emits per condition, into cn_bar/:
#   <exp>_<cond>_cn_bar.csv    Stage 4 contract: CHROM, POS, cn_bar
#   <exp>_<cond>_audit.tsv     per-SV status and the breakpoint cn_bars used
#
# Usage: ./09_chrx_join_all.sh [--weighted]

set -euo pipefail

_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRA=("$@")

# For cov and the coverage-named truth/conditions paths. Sourced rather than
# re-deriving "${COV:-25}" here, so this cannot drift from the script that writes
# them.  Run as:  COV=50 SEG_DIR=$V/seg10k OUT_DIR=$V/cn_bar10k_m ./09_...sh
source "${_DIR}/chrx_common.sh"

# SEG_DIR / OUT_DIR are overridable so a join against a finer re-segmentation
# does not clobber the first one:
#   SEG_DIR=$V/seg10k OUT_DIR=$V/cn_bar10k ./09_chrx_join_all.sh
truth="$truth_tsv"
conds="$cond_tsv"
seg_dir="${SEG_DIR:-${BASE_DIR}/seg}${rep_tag}"      # rep-safe; see chrx_common.sh
out_dir="${OUT_DIR:-${BASE_DIR}/cn_bar}${rep_tag}"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$out_dir"

[[ -s "$truth" ]] || { echo "ERROR: missing $truth (run 04_chrx_ground_truth.sh)" >&2; exit 1; }

n=0
while IFS=$'\t' read -r _tid exp cond _rest; do
    seg="${seg_dir}/chrx_seg_${exp}_${cond}.csv"
    if [[ ! -s "$seg" ]]; then
        echo "SKIP ${exp}/${cond}: no segmentation at $seg" >&2
        continue
    fi
    awk -F'\t' -v e="$exp" -v c="$cond" 'NR==1 || ($1==e && $2==c)' "$truth" > "$tmp/svs.tsv"
    python3 "${_DIR}/08_chrx_segment_sv_join.py" \
        --segments "$seg" --svs "$tmp/svs.tsv" \
        --out   "${out_dir}/${exp}_${cond}_cn_bar.csv" \
        --audit "${out_dir}/${exp}_${cond}_audit.tsv" \
        ${EXTRA[@]+"${EXTRA[@]}"}
    n=$((n + 1))
done < <(tail -n +2 "$conds")

echo
echo "RESULT: joined $n condition(s) -> $out_dir"
