#!/bin/bash
#SBATCH --job-name=vis_mhack
#SBATCH --output=log/hack_%A_%a.out
#SBATCH --error=log/hack_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --time=2:00:00
#SBATCH --array=0


BASE_DIR="${BASE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ref="${REF:-/path/to/hg38_chr1-2.fa}"
in_dir="${BASE_DIR}/resources/beds"
sv_beds_dir="${BASE_DIR}/truth/sv_beds"
snp_base="${BASE_DIR}/truth/fastas/snp_background"
fasta_out="${BASE_DIR}/truth/fastas/clone_genomes"
scrp_dir="${BASE_DIR}/scripts/genome_setup"
clone_lst=(c2 c3)
cc_lst=(c22 c33)
cln=${clone_lst[$SLURM_ARRAY_TASK_ID]}
cc=${cc_lst[$SLURM_ARRAY_TASK_ID]}

if [ -f "$sv_beds_dir/c22.bed" ]; then
    echo "File exists. Skipping command."
else
    echo "File does not exist. Running command..."
    Rscript $scrp_dir/sv2subclone.R \
        -a $in_dir/H1.bed \
        -b $in_dir/H1_1.bed \
        -o $sv_beds_dir
fi



# the first step where SNP is in genome, use data from single clone
## exp1: SV alone
#VISOR HACk -g $ref -b $in_dir/chr1_2_snps.bed -o $snp_base/snp1 #put het snp on
VISOR HACk -g $snp_base/snp1/h1.fa -b $sv_beds_dir/$cln.bed -o $fasta_out/${cln}a1 #put het sv on
VISOR HACk -g $ref -b $sv_beds_dir/hom_$cln.bed -o $fasta_out/${cln}a2 #put hom sv on
mkdir -p $fasta_out/e1/$cln
mv $fasta_out/${cln}a1/* $fasta_out/e1/$cln/
mv $fasta_out/${cln}a2/h1.fa $fasta_out/e1/$cln/h2.fa
mv $fasta_out/${cln}a2/h1.fa.fai $fasta_out/e1/$cln/h2.fa.fai
rm -r $fasta_out/${cln}a1 $fasta_out/${cln}a2

## exp2: sv before in-cis amplification
#VISOR HACk -g $ref -b $in_dir/o_chr1_2_snps.bed -o $snp_base/snp1_1
VISOR HACk -g $snp_base/snp1_1/h1.fa -b $sv_beds_dir/$cc.bed -o $fasta_out/tmp_$cln #put het sv on
VISOR HACk -g $ref -b $sv_beds_dir/hom_$cc.bed -o $fasta_out/${cln}a2 #put hom sv on
VISOR HACk -g $fasta_out/tmp_$cln/h1.fa -b $sv_beds_dir/dup_$cc.bed -o $fasta_out/${cln}a1 #put in-cis amp
mkdir -p $fasta_out/e2/$cln
mv $fasta_out/${cln}a1/* $fasta_out/e2/$cln/
mv $fasta_out/${cln}a2/h1.fa $fasta_out/e2/$cln/h2.fa
mv $fasta_out/${cln}a2/h1.fa.fai $fasta_out/e2/$cln/h2.fa.fai
rm -r $fasta_out/${cln}a1 $fasta_out/${cln}a2

## exp3: sv before in-trans amplification
# use existing tmp from exp2 (has snp and het SV) as ale1
##clone1
VISOR HACk -g $ref -b $sv_beds_dir/hom_$cc.bed -o $fasta_out/itmd #put hom sv on
VISOR HACk -g $fasta_out/itmd/h1.fa -b $sv_beds_dir/dup_$cc.bed -o $fasta_out/${cln}a2 #put in-trans amp
mkdir -p $fasta_out/e3/$cln
cp $fasta_out/tmp_$cln/* $fasta_out/e3/$cln
mv $fasta_out/${cln}a2/h1.fa $fasta_out/e3/$cln/h2.fa
mv $fasta_out/${cln}a2/h1.fa.fai $fasta_out/e3/$cln/h2.fa.fai
rm -r $fasta_out/itmd $fasta_out/${cln}a2


## exp4: sv after in-cis amplification
# use existing snp1_1 to start (with het snp)
##clone1
VISOR HACk -g $snp_base/snp1_1/h1.fa -b $sv_beds_dir/dup_$cc.bed -o $fasta_out/itmd #put amp on
VISOR HACk -g $fasta_out/itmd/h1.fa -b $sv_beds_dir/aft_$cc.bed -o $fasta_out/${cln}a1 # put het sv on
VISOR HACk -g $ref -b $sv_beds_dir/hom_$cc.bed -o $fasta_out/${cln}a2 # put hom sv on
mkdir -p $fasta_out/e4/$cln
mv $fasta_out/${cln}a1/* $fasta_out/e4/$cln
mv $fasta_out/${cln}a2/h1.fa $fasta_out/e4/$cln/h2.fa
mv $fasta_out/${cln}a2/h1.fa.fai $fasta_out/e4/$cln/h2.fa.fai
rm -r $fasta_out/itmd $fasta_out/${cln}a2 $fasta_out/${cln}a1

## exp5: sv after del
# use tmp for allele1
##clone1,2,3
VISOR HACk -g $ref -b $sv_beds_dir/del_$cc.bed -o $fasta_out/a2
mkdir -p $fasta_out/e5/$cln
cp $fasta_out/tmp_$cln/* $fasta_out/e5/$cln
mv $fasta_out/a2/h1.fa $fasta_out/e5/$cln/h2.fa
mv $fasta_out/a2/h1.fa.fai $fasta_out/e5/$cln/h2.fa.fai
rm -r $fasta_out/a2

cp $ref $snp_base/snp1/h2.fa
cp $ref.fai $snp_base/snp1/h2.fa.fai
cp $ref $snp_base/snp1_1/h2.fa
cp $ref.fai $snp_base/snp1_1/h2.fa.fai


