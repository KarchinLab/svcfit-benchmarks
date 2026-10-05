#!/bin/bash
# Prepare chromosome 1 and 2 SNP VCFs from dbSNP build 157 (see download_snp.sh):
# rename chromosome 1 RefSeq accessions to chr names and extract chromosome 2.
bcftools annotate --rename-chrs chrom_map.txt -Oz -o chr1_renamed.vcf.gz chr1.vcf.gz
bcftools view -r NC_000002.12 GCF_000001405.40.gz -Oz -o chr2.vcf.gz

