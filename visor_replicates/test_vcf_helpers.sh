#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
rscript=${R_SCRIPT:-$(command -v Rscript || true)}
[[ -n "$rscript" && -x "$rscript" ]] || {
    echo 'ERROR: set R_SCRIPT to an Rscript with optparse, readr, and dplyr' >&2
    exit 1
}

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/visor-vcf-helpers.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT
mkdir -p "$tmp_dir/ranges" "$tmp_dir/svclone" "$tmp_dir/facet"

vcf="$tmp_dir/input.vcf"
printf '%s\n' \
    '##fileformat=VCFv4.2' \
    $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tnormal\ttumor' \
    $'chr1\t100\tbnd1\tN\tN]chr2:200]\t60\tPASS\tSVTYPE=BND\tGT\t0/0\t0/1' \
    $'chr1\t1000\tdel1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1200\tGT\t0/0\t0/1' \
    $'chr2\t2000\tdup1\tN\t<DUP:TANDEM>\t60\tPASS\tSVTYPE=DUP;END=2400\tGT\t0/0\t0/1' \
    $'chr2\t3000\tinv1\tN\t<INV>\t60\tPASS\tSVTYPE=INV;END=3500\tGT\t0/0\t0/1' \
    $'chr1\t4000\tins1\tN\t<INS>\t60\tPASS\tSVTYPE=INS;END=4001\tGT\t0/0\t0/1' \
    > "$vcf"

"$rscript" "$script_dir/helpers/get_sv_range.R" \
    -s fixture -o "$tmp_dir/ranges" -v "$vcf"
ranges="$tmp_dir/ranges/fixture.bed"
test "$(wc -l < "$ranges")" -eq 10
awk -F '\t' 'NF != 3 || $1 == "" || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ { bad=1 } END { exit bad }' "$ranges"

TASK_FIXTURE_ROOT="$tmp_dir" "$rscript" -e '
root <- Sys.getenv("TASK_FIXTURE_ROOT")
fit <- list(ploidy = 2)
save(fit, file = file.path(root, "facet", "fixture.RData"))
segments <- data.frame(
  chrom = c(1, 2), start = c(1, 1), end = c(5000000, 5000000),
  tcn.em = c(2, 2), lcn.em = c(1, 1)
)
write.table(
  segments, file.path(root, "facet", "fixture.bed"), sep = "\t",
  row.names = FALSE, quote = FALSE
)
'

"$rscript" "$script_dir/helpers/make_input.R" \
    -p 0.1 -f "$vcf" -s fixture -o "$tmp_dir/svclone" -c "$tmp_dir/facet"
simple="$tmp_dir/svclone/svc_fixture_simple.txt"
test "$(wc -l < "$simple")" -eq 5
grep -Fxq $'chr1\tpos1\tchr2\tpos2' "$simple"
awk -F '\t' 'NR > 1 && (NF != 4 || $1 == "" || $2 !~ /^[0-9]+$/ || $3 == "" || $4 !~ /^[0-9]+$/) { bad=1 } END { exit bad }' "$simple"
grep -Fxq $'chr1\t100\tchr2\t200' "$simple"
grep -Fxq $'chr1\t1000\tchr1\t1200' "$simple"
grep -Fxq $'chr2\t2000\tchr2\t2400' "$simple"
grep -Fxq $'chr2\t3500\tchr2\t3000' "$simple"
test -s "$tmp_dir/svclone/spp_fixture.txt"
test -s "$tmp_dir/svclone/cnv.txt"

sed 's/SVTYPE=DEL;END=1200/SVTYPE=DEL/' "$vcf" > "$tmp_dir/missing-end.vcf"
if "$rscript" "$script_dir/helpers/make_input.R" \
    -p 0.1 -f "$tmp_dir/missing-end.vcf" -s fixture \
    -o "$tmp_dir/invalid" -c "$tmp_dir/facet" >/dev/null 2>&1; then
    echo 'ERROR: retained DEL without END unexpectedly passed' >&2
    exit 1
fi

grep -Fxq 'mean_cov: 50' "$script_dir/helpers/svclone_config.ini"
grep -Fxq 'chroms: 1,2' "$script_dir/helpers/svclone_config.ini"
grep -Fxq 'clus_limit: 3' "$script_dir/helpers/svclone_config.ini"

echo 'PASS: VISOR VCF helper compatibility'
