#!/bin/bash

BASE_DIR="${BASE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
# REF_DIR: directory containing individual chromosome FASTAs (chr1.fa, chr2.fa)
REF_DIR="${REF_DIR:-/path/to/hg38_reference/indv_chr}"
sv_size=100000
sv_num=110
use_sv=100
input_dir="${BASE_DIR}/resources/beds"
scrp_dir="${BASE_DIR}/scripts/genome_setup"
rm $input_dir/*.bed


## sv
curl -L -o $input_dir/sv_repeat_telomere_centromere.bed https://gist.githubusercontent.com/chapmanb/4c40f961b3ac0a4a22fd/raw/2025f3912a477edc597e61d911bd1044dc943440/sv_repeat_telomere_centromere.bed

echo -e "chr1" > $input_dir/tmp.txt && grep -w -f $input_dir/tmp.txt $input_dir/sv_repeat_telomere_centromere.bed \
	| sortBed > $input_dir/exclude1.bed && rm $input_dir/tmp.txt
echo -e "chr2" > $input_dir/tmp.txt && grep -w -f $input_dir/tmp.txt $input_dir/sv_repeat_telomere_centromere.bed \
	| sortBed > $input_dir/exclude2.bed && rm $input_dir/tmp.txt

cut -f1,2 ${REF_DIR}/chr1.fa.fai > $input_dir/chrom1.dim.tsv
cut -f1,2 ${REF_DIR}/chr2.fa.fai > $input_dir/chrom2.dim.tsv

## for exp1 where no overlaps
Rscript ${VISOR_HOME:?set VISOR_HOME to the VISOR checkout}/scripts/randomregion.r\
        -d $input_dir/chrom1.dim.tsv\
        -n $sv_num\
        -l $sv_size\
        -x $input_dir/exclude1.bed\
        -v 'deletion,inversion,tandem duplication,reciprocal translocation'\
        -r '25:25:25:25'\
        | sortBed > $input_dir/HACk1.bed

## for sv overlaps with cnv (hence no duplication)
Rscript ${VISOR_HOME}/scripts/randomregion.r\
        -d $input_dir/chrom1.dim.tsv\
        -n $sv_num\
        -l $sv_size\
        -x $input_dir/exclude1.bed\
        -v 'deletion,inversion,reciprocal translocation'\
        -r '34:33:33'\
        | sortBed > $input_dir/HACk11.bed

## select locations for reciprocal trans on chr2
Rscript ${VISOR_HOME}/scripts/randomregion.r\
        -d $input_dir/chrom2.dim.tsv\
        -n $sv_num\
        -l $sv_size\
        -x $input_dir/exclude2.bed\
        -v 'deletion'\
        -r '100'\
        | sortBed > $input_dir/H2.bed

## create inter-chrom reciprocal trans using chr1 and chr2
Rscript $scrp_dir/modify_sv.R \
	-a $input_dir/HACk1.bed \
	-b $input_dir/H2.bed \
	-n $use_sv \
	-z $sv_size \
	-o $input_dir

Rscript $scrp_dir/modify_sv.R \
	-a $input_dir/HACk11.bed \
	-b $input_dir/H2.bed \
	-n $use_sv \
	-z $sv_size \
	-o $input_dir

rm $input_dir/exclude* $input_dir/*.tsv $input_dir/sv_repeat_telomere_centromere.bed

## select het SNP based on SV
bash $scrp_dir/filter_snp.sh "${BASE_DIR}/resources/snp_vcfs"

## create overlap cnv (also create aH1_1.bed, adjust sv location based on dup happened before SV)
Rscript $scrp_dir/make_cnv.R -n $use_sv -p $input_dir -o $input_dir

#total=$(wc -l < H1_1.bed)
#half=$(( (total + 1) / 2 ))   # ceil(total/2)
#head -n "$half" H1_1.bed > del_h11.bed
#tail -n +$((half + 1)) H1_1.bed > del2_h11.bed
