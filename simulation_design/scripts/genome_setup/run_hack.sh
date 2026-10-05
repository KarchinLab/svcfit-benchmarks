#!/bin/bash
#SBATCH --job-name=vis_hack
#SBATCH --output=log/hack_%A_%a.out
#SBATCH --error=log/hack_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=10:00:00

BASE_DIR="${BASE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ref="${REF:-/path/to/hg38_chr1-2.fa}"
in_dir="${BASE_DIR}/resources/beds"
snp_base="${BASE_DIR}/truth/fastas/snp_background"
fasta_out="${BASE_DIR}/truth/fastas/clone_genomes"

#bash $scrp_dir/make_sv.sh 

##exp1: sv alone
echo --------starting e1--------
#VISOR HACk -g $ref -b $in_dir/chr1_2_snps.bed -o $fasta_out/snp1 #put het snp on
VISOR HACk -g $snp_base/snp1/h1.fa -b $in_dir/H1.bed -o $fasta_out/ale1 #put het sv on
VISOR HACk -g $ref -b $in_dir/hom_H1.bed -o $fasta_out/ale2 #put hom sv on
mkdir -p $fasta_out/e1
mv $fasta_out/ale1/* $fasta_out/e1
mv $fasta_out/ale2/h1.fa $fasta_out/e1/h2.fa
mv $fasta_out/ale2/h1.fa.fai $fasta_out/e1/h2.fa.fai
rm -r $fasta_out/ale1 $fasta_out/ale2


## exp2: sv before in-cis amplification
echo --------starting e2--------
#VISOR HACk -g $ref -b $in_dir/o_chr1_2_snps.bed -o $fasta_out/snp1_1 #put het snp on
VISOR HACk -g $snp_base/snp1_1/h1.fa -b $in_dir/H1_1.bed -o $fasta_out/tmp1 #put het sv on
VISOR HACk -g $ref -b $in_dir/hom_H1_1.bed -o $fasta_out/ale2 #put hom sv on
VISOR HACk -g $fasta_out/tmp1/h1.fa -b $in_dir/dup_over.bed -o $fasta_out/ale1 #put in-cis amp
mkdir -p $fasta_out/e2
mv $fasta_out/ale1/* $fasta_out/e2
mv $fasta_out/ale2/h1.fa $fasta_out/e2/h2.fa
mv $fasta_out/ale2/h1.fa.fai $fasta_out/e2/h2.fa.fai
rm -r $fasta_out/ale1 $fasta_out/ale2

## exp3: sv before in-trans amplification
echo --------starting e3--------
# use existing tmp1 from exp2 (has snp and het SV) as ale1
VISOR HACk -g $ref -b $in_dir/hom_H1_1.bed -o $fasta_out/tmp2 #put hom sv on
VISOR HACk -g $fasta_out/tmp2/h1.fa -b $in_dir/dup_over.bed -o $fasta_out/ale2 #put in-cis amp
mkdir -p $fasta_out/e3
cp $fasta_out/tmp1/* $fasta_out/e3
mv $fasta_out/ale2/h1.fa $fasta_out/e3/h2.fa
mv $fasta_out/ale2/h1.fa.fai $fasta_out/e3/h2.fa.fai
rm -r $fasta_out/tmp2 $fasta_out/ale2

## exp4: sv after in-cis amplification
echo --------starting e4--------
# use existing snp1_1 to start (with het snp)
VISOR HACk -g $snp_base/snp1_1/h1.fa -b $in_dir/dup_over.bed -o $fasta_out/tmp11 #put amp on
VISOR HACk -g $fasta_out/tmp11/h1.fa -b $in_dir/aH1_1.bed -o $fasta_out/ale1 # put het sv on
VISOR HACk -g $ref -b $in_dir/hom_H1_1.bed -o $fasta_out/ale2 # put hom sv on
mkdir -p $fasta_out/e4
cp $fasta_out/ale1/* $fasta_out/e4
mv $fasta_out/ale2/h1.fa $fasta_out/e4/h2.fa
mv $fasta_out/ale2/h1.fa.fai $fasta_out/e4/h2.fa.fai
rm -r $fasta_out/ale2 $fasta_out/tmp11 $fasta_out/ale1


## exp5: sv after del
echo --------starting e5--------
VISOR HACk -g $ref -b $in_dir/del1.bed -o $fasta_out/ale2
mkdir -p $fasta_out/e5
cp $fasta_out/tmp1/* $fasta_out/e5
mv $fasta_out/ale2/h1.fa $fasta_out/e5/h2.fa
mv $fasta_out/ale2/h1.fa.fai $fasta_out/e5/h2.fa.fai
rm -r $fasta_out/tmp1 $fasta_out/ale2

## mv ref to SNP to create a ref with het snp
#cp $ref $fasta_out/snp1/h2.fa
#cp $ref.fai $fasta_out/snp1/h2.fai
#cp $ref $fasta_out/snp1_1/h2.fa
#cp $ref.fai $fasta_out/snp1_1/h2.fai
