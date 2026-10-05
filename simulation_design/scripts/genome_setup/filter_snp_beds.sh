#!/bin/bash
# One-time script: filter SNP BEDs to ±1000 bp around SV breakpoints.
# Outputs are written to the same directory and are read by the calling stages.

res_dir="${1:?usage: filter_snp_beds.sh <resources/beds directory>}"
bedtools="${BEDTOOLS:-bedtools}"

make_windows() {
    local sv_beds=("$@")
    awk 'BEGIN{OFS="\t"} !/^#/ {
        s=$2-1000; if(s<0) s=0; print $1, s, $2+1000
        e=$3-1000; if(e<0) e=0; print $1, e, $3+1000
    }' "${sv_beds[@]}" | sort -k1,1 -k2,2n | $bedtools merge -i stdin
}

# exp1: H1.bed (chr1) + H2.bed (chr2)
make_windows $res_dir/H1.bed $res_dir/H2.bed \
    | $bedtools intersect -a $res_dir/chr1_2_snps.bed -b stdin \
    > $res_dir/chr1_2_snps_sv1kb.bed
echo "Written: chr1_2_snps_sv1kb.bed ($(wc -l < $res_dir/chr1_2_snps_sv1kb.bed) SNPs)"

# exp2-5: H1_1.bed (chr1) + H2.bed (chr2)
make_windows $res_dir/H1_1.bed $res_dir/H2.bed \
    | $bedtools intersect -a $res_dir/o_chr1_2_snps.bed -b stdin \
    > $res_dir/o_chr1_2_snps_sv1kb.bed
echo "Written: o_chr1_2_snps_sv1kb.bed ($(wc -l < $res_dir/o_chr1_2_snps_sv1kb.bed) SNPs)"
