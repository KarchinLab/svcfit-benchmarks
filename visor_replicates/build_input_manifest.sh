#!/usr/bin/env bash

# Build the complete immutable-input mapping for a 2,250-condition VISOR run.
# The output is provenance and an executable interface: workers verify their
# decoded source/result paths against the corresponding manifest row.

set -euo pipefail
umask 0027

source_root=""
result_root=""
output=""

usage() {
    echo "Usage: $0 --source-root DIR --result-root DIR --output FILE"
}

while (( $# )); do
    case "$1" in
        --source-root) source_root=${2:?}; shift 2 ;;
        --result-root) result_root=${2:?}; shift 2 ;;
        --output) output=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

source_root=${source_root%/}
result_root=${result_root%/}
[[ "$source_root" = /* && -d "$source_root" ]] || {
    echo "ERROR: --source-root must be an existing absolute directory" >&2; exit 2;
}
[[ "$result_root" = /* && "$result_root" != "$source_root" ]] || {
    echo "ERROR: --result-root must be absolute and differ from --source-root" >&2; exit 2;
}
[[ "$output" = /* ]] || {
    echo "ERROR: --output must be an absolute path" >&2; exit 2;
}
[[ ! -L "$output" ]] || {
    echo "ERROR: manifest path must not be a symbolic link: $output" >&2; exit 3;
}
[[ ! -e "$output" || -f "$output" ]] || {
    echo "ERROR: manifest path exists and is not a regular file: $output" >&2; exit 3;
}
mkdir -p "$(dirname "$output")"

tmp="${output}.inprogress.$$"
printf 'task_id\treplicate\texperiment\tcondition\tsource_work\tresult_work\tshort_dir\tsnp_dir\tfacet_dir\tmanta_vcf\tsvtyp_vcf\n' > "$tmp"

purities=(10 10 10 20 20 20 40 40 40 60 60 60 80 80 80)
mixtures=(10 30 50 10 30 50 10 30 50 10 30 50 10 30 50)
for task_id in $(seq 0 2249); do
    rep_idx=$(( task_id / 75 ))
    remainder=$(( task_id % 75 ))
    exp_idx=$(( remainder / 15 ))
    cond_idx=$(( remainder % 15 ))
    rep=rep$((rep_idx + 1))
    experiment=exp$((exp_idx + 1))
    condition=c50p${purities[$cond_idx]}m${mixtures[$cond_idx]}
    source_work="$source_root/$rep/$experiment/$condition"
    result_work="$result_root/$rep/$experiment/$condition"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$task_id" "$rep" "$experiment" "$condition" \
        "$source_work" "$result_work" \
        "$source_work/short" "$source_work/SNP" "$source_work/facet" \
        "$source_work/manta/results/variants/${condition}.vcf" \
        "$source_work/svtyp/svt_${condition}.vcf" >> "$tmp"
done

[[ "$(wc -l < "$tmp")" == 2251 ]] || {
    echo "ERROR: manifest does not contain header plus 2,250 rows" >&2
    exit 4
}
if [[ -e "$output" ]]; then
    if ! cmp -s "$tmp" "$output"; then
        rm -f -- "$tmp"
        echo "ERROR: existing manifest does not match requested roots: $output" >&2
        exit 3
    fi
    rm -f -- "$tmp"
    echo "Verified existing immutable-input manifest: $output"
else
    mv "$tmp" "$output"
    echo "Wrote immutable-input manifest: $output"
fi
