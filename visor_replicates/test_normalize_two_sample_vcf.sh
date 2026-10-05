#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
normalizer="$script_dir/normalize_two_sample_vcf.awk"
tmp_dir=$(mktemp -d "$script_dir/.normalize-vcf-test.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT

input="$tmp_dir/input.vcf"
output="$tmp_dir/output.vcf"

printf '%s\n' \
    '##fileformat=VCFv4.2' \
    '##source=test' \
    $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbulk\tbulk' \
    $'chr1\t10\tbnd1\tN\tN]chr2:20]\t60\tPASS\tSVTYPE=BND\tGT\t0/0\t0/1' \
    > "$input"

awk \
    -v normal_sample=rep1_exp1_c50p10m10_normal \
    -v tumor_sample=rep1_exp1_c50p10m10_tumor \
    -f "$normalizer" "$input" > "$output"

grep -Fxq '##SVCFitSampleNameNormalization=<Normal=rep1_exp1_c50p10m10_normal,Tumor=rep1_exp1_c50p10m10_tumor>' "$output"
grep -Fxq $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\trep1_exp1_c50p10m10_normal\trep1_exp1_c50p10m10_tumor' "$output"
test "$(awk '!/^#/ { n++ } END { print n + 0 }' "$output")" -eq 1
grep -Fxq $'chr1\t10\tbnd1\tN\tN]chr2:20]\t60\tPASS\tSVTYPE=BND\tGT\t0/0\t0/1' "$output"

# Exercise the temporary BAM-name compatibility round trip used around
# SVtyper. It must retain exactly one metadata line and end with two unique,
# deterministic samples.
awk \
    -v normal_sample=rep1_exp1_c50p10m10_normal \
    -v tumor_sample=bulk \
    -v emit_metadata=0 \
    -f "$normalizer" "$output" > "$tmp_dir/compat.vcf"
grep -Fxq $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\trep1_exp1_c50p10m10_normal\tbulk' "$tmp_dir/compat.vcf"

awk \
    -v normal_sample=rep1_exp1_c50p10m10_normal \
    -v tumor_sample=rep1_exp1_c50p10m10_tumor \
    -v emit_metadata=0 \
    -f "$normalizer" "$tmp_dir/compat.vcf" > "$tmp_dir/final.vcf"
grep -Fxq $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\trep1_exp1_c50p10m10_normal\trep1_exp1_c50p10m10_tumor' "$tmp_dir/final.vcf"
test "$(grep -Fc '##SVCFitSampleNameNormalization=' "$tmp_dir/final.vcf")" -eq 1
test "$(awk -F '\t' '!/^#/ { print NF; exit }' "$tmp_dir/final.vcf")" -eq 11

printf '%s\n' \
    '##fileformat=VCFv4.2' \
    $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tonly_one' \
    > "$tmp_dir/one-sample.vcf"
if awk -v normal_sample=n -v tumor_sample=t -f "$normalizer" "$tmp_dir/one-sample.vcf" >/dev/null 2>&1; then
    echo 'ERROR: one-sample VCF unexpectedly passed' >&2
    exit 1
fi

if awk -v normal_sample=same -v tumor_sample=same -f "$normalizer" "$input" >/dev/null 2>&1; then
    echo 'ERROR: duplicate replacement names unexpectedly passed' >&2
    exit 1
fi

printf '%s\n' \
    '##fileformat=VCFv4.2' \
    $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tnormal\ttumor\tbulk' \
    $'chr1\t10\tbnd1\tN\tN]chr2:20]\t60\tPASS\tSVTYPE=BND\tGT\t0/0\t0/1\t0/1' \
    > "$tmp_dir/three-sample.vcf"
if awk -v normal_sample=n -v tumor_sample=t -f "$normalizer" "$tmp_dir/three-sample.vcf" >/dev/null 2>&1; then
    echo 'ERROR: three-sample VCF unexpectedly passed' >&2
    exit 1
fi

printf '%s\n' '##fileformat=VCFv4.2' > "$tmp_dir/no-header.vcf"
if awk -v normal_sample=n -v tumor_sample=t -f "$normalizer" "$tmp_dir/no-header.vcf" >/dev/null 2>&1; then
    echo 'ERROR: headerless VCF unexpectedly passed' >&2
    exit 1
fi

echo 'PASS: normalize_two_sample_vcf.awk'
