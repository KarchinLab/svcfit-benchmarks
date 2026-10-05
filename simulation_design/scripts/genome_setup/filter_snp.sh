
# 1. Ensure exactly one argument was provided
if [ $# -ne 1 ]; then
  echo "Usage: $0 <pathway>"
  exit 1
fi

# 2. Capture the argument into a descriptive variable
dir="$1"

if [ ! -e "$dir" ]; then
  echo "Error: Path '$dir' does not exist."
  exit 2
fi

## filter snps to be within SV breakpoints
#1. generate sv interval for non-overlap
awk 'BEGIN{OFS="\t"}
{
  chr  = $1
  s    = $2     # SV start
  e    = $3     # SV end
  # window around start
  ws   = (s-30000>0 ? s-30000 : 0)
  we   = s+30000
  print chr, ws, we
  # window around end
  es   = (e-30000>0 ? e-30000 : 0)
  ee   = e+30000
  print chr, es, ee
}' $dir/../beds/H1.bed > $dir/sv_flank1.bed

bcftools view \
  -R $dir/sv_flank1.bed \
  -Oz \
  -o  $dir/chr1_svsnp.vcf.gz \
  $dir/chr1_snps.vcf.gz

tabix -p vcf $dir/chr1_svsnp.vcf.gz

## for overlaps
awk 'BEGIN{OFS="\t"}
{
  chr  = $1
  s    = $2     # SV start
  e    = $3     # SV end
  # window around start
  ws   = (s-30000>0 ? s-30000 : 0)
  we   = s+30000
  print chr, ws, we
  # window around end
  es   = (e-30000>0 ? e-30000 : 0)
  ee   = e+30000
  print chr, es, ee
}' $dir/../beds/H1_1.bed > $dir/sv_flank2.bed

bcftools view \
  -R $dir/sv_flank2.bed \
  -Oz \
  -o  $dir/chr1_1_svsnp.vcf.gz \
  $dir/chr1_snps.vcf.gz

tabix -p vcf $dir/chr1_1_svsnp.vcf.gz


## for chr2
awk 'BEGIN{OFS="\t"}
{
  chr  = $1
  s    = $2     # SV start
  e    = $3     # SV end
  # window around start
  ws   = (s-30000>0 ? s-30000 : 0)
  we   = s+30000
  print chr, ws, we
  # window around end
  es   = (e-30000>0 ? e-30000 : 0)
  ee   = e+30000
  print chr, es, ee
}' $dir/../beds/H2.bed > $dir/sv_flank3.bed

bcftools view \
  -R $dir/sv_flank3.bed \
  -Oz \
  -o  $dir/chr2_svsnp.vcf.gz \
  $dir/chr2_snps.vcf.gz

tabix -p vcf $dir/chr2_svsnp.vcf.gz


## formatting
bcftools query   -f '%CHROM\t%POS\t%REF\t%ALT\n'   $dir/chr1_svsnp.vcf.gz > $dir/chr1_snps.txt

bcftools query   -f '%CHROM\t%POS\t%REF\t%ALT\n'   $dir/chr1_1_svsnp.vcf.gz > $dir/chr1_1_snps.txt

bcftools query   -f '%CHROM\t%POS\t%REF\t%ALT\n'   $dir/chr2_svsnp.vcf.gz > $dir/chr2_snps.txt

## get file into VISOR input format
awk -F'\t' 'BEGIN{OFS="\t"}
{
	split($4, a, ",")
	base = substr(a[1], 1, 1)
	print $1, $2-1, $2, "SNP", base, 0
}' $dir/chr1_snps.txt > $dir/chr1_snps.full

awk -F'\t' 'BEGIN{OFS="\t"}
{
	split($4, a, ",")
	base = substr(a[1], 1, 1)
	print $1, $2-1, $2, "SNP", base, 0
}' $dir/chr1_1_snps.txt > $dir/chr1_1_snps.full

awk -F'\t' 'BEGIN{OFS="\t"}
{
        split($4, a, ",")
        base = substr(a[1], 1, 1)
        print $1, $2-1, $2, "SNP", base, 0
}' $dir/chr2_snps.txt > $dir/chr2_snps.full

awk 'NR % 30 == 1' $dir/chr1_snps.full > $dir/tmp_chr1_snps.bed
awk 'NR % 30 == 1' $dir/chr1_1_snps.full > $dir/tmp_chr1_1_snps.bed
awk 'NR % 30 == 1' $dir/chr2_snps.full > $dir/tmp_chr2_snps.bed

cat $dir/tmp_chr1_snps.bed $dir/tmp_chr2_snps.bed > $dir/../beds/chr1_2_snps.bed
cat $dir/tmp_chr1_1_snps.bed $dir/tmp_chr2_snps.bed > $dir/../beds/o_chr1_2_snps.bed


rm $dir/*.txt $dir/*.full $dir/chr1_svsnp.vcf.gz $dir/sv_flank2.bed $dir/sv_flank1.bed $dir/sv_flank3.bed
