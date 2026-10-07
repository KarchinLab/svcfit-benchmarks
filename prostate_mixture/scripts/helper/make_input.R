library(dplyr)
library(optparse)

option_list <- list(
  make_option(c("-f", "--file"),
              type = "character",
              help = "SV vcf file"),
  #make_option(c("-v", "--cov"),
  #            type = "character",
  #            help = "cov number"),
  make_option(c("-p", "--purity"),
              type = "character",
              help = "bam purity"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "output dir"),
  #make_option(c("-k", "--hack"),
  #            type = "character",
  #            help = "hack dir"),
  make_option(c("-c", "--cnv"),
              type = "character",
              help = "facet dir"),
  make_option(c("-s", "--smp"),
              type = "character",
              help = "sample name")
)

parser <- OptionParser(option_list = option_list,add_help_option = FALSE)
opts   <- parse_args(parser)

#samp=paste0("c",opts$cov,"p",opts$purity)
samp=opts$smp
path=opts$file
svt <- read.table(path, quote="\"")
colnames(svt) <- c("CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO",    "FORMAT","normal","tumor")
tmp <- svt %>%
      dplyr::rename("chr1" = "CHROM", "pos1" = "POS")%>%
      filter(FILTER == "PASS")%>%
      ## BND mate chromosome, fixed 2026-10-06. The old regex ".*(chr.).*" only matched a
      ## "chr"-prefixed mate, but these hs37d5 Manta VCFs write mates as e.g. "[7:45702526[",
      ## so every translocation fell back to chr2 = chr1 and became a fake intrachromosomal
      ## event. Take the mate contig from the bracket notation, whatever its naming.
      mutate(chr2 = ifelse(grepl("[][]", ALT),
                           sub("^[^][]*[][]([^:]+):.*$", "\\1", ALT), chr1),
             pos2 = gsub(".*:(\\d+).*","\\1", ALT),
             ## END= anchored to a key boundary (2026-10-06): SVtyper rewrites INFO so it no longer
             ## starts with END=, which the old unanchored gsub turned into NA (or matched CIEND=).
             pos2 = ifelse(grepl("\\d+", pos2), pos2, sub("^(.*;)?END=(\\d+).*$","\\2",INFO)),
             pos2 = as.integer(pos2),
             class = gsub(".*SVTYPE=(\\w+).*","\\1", INFO),
             pp=pos1,
             pos1=ifelse(class=="INV", pos2, pos1),
             pos2=ifelse(class=="INV", pp, pos2))%>%
      filter(!class=="INS") %>%
      select(chr1,pos1,chr2,pos2)

write.table(tmp, paste0(opts$output,"/","svc_",samp,"_simple.txt"),row.names = F,col.names = T, sep="\t", quote=FALSE)

load(paste0(opts$cnv, "/", samp,".RData"))
ploi=fit$ploidy
pur=ifelse(is.na(fit$purity), as.numeric(opts$purity), fit$purity)
pp = data.frame(sample=samp, purity=pur, ploidy = ploi)
write.table(pp, paste0(opts$output,"/","spp_",samp,".txt"),row.names = F,col.names = T, sep="\t", quote=FALSE)

###################### make cnv #########################
facet <- read.delim(paste0(opts$cnv,"/",samp,".bed"))
dat = facet %>%
  mutate(chrom = as.character(chrom),
         ## chrX naming. FACETS codes chrX as 23 (24 for Y), but the SVs -- from Manta via
         ## simple.txt -- use "X". SVclone matches an SV to its copy-number segment by chromosome
         ## name, so emit "X"/"Y" so the CNV matches.
         chrom = ifelse(chrom == "23", "X", ifelse(chrom == "24", "Y", chrom)),
         n_major=2,
         n_minor=1,
         t_total=tcn.em,
         t_minor=ifelse(is.na(lcn.em), 0L, lcn.em))%>%
  select(chrom, start, end, n_major, n_minor, t_total, t_minor)
#SV_del=read.delim(paste0(opts$hack,"/SV_del.bed"), header=FALSE)
#SV_dup=read.delim(paste0(opts$hack,"/SV_dup.bed"), header=FALSE)
#SV_cnv=rbind(SV_dup, SV_del)
#SV_dup2=read.delim(paste0(opts$hack,"/SV_dup2.bed"), header=FALSE)
#colnames(SV_cnv) <- c("CHROM","start","end","type","info","flank")
#colnames(SV_dup2) <- c("CHROM","start","end","type","info","flank")
#cnv=rbind(SV_cnv,SV_dup2) %>%
#  mutate(norm_tcn=2,
#         norm_minor=1,
#         tum_tcn=ifelse(type=="deletion", 1, 3),
#         tum_minor=ifelse(type=="deletion", 0, 1),
#         CHROM="chr22")%>%
#  select(-c(type, info, flank))%>%
#  rbind(c("chr22",1,(min(SV_cnv$start)-1),2,1,2,1))%>%
#  arrange(start)

write.table(dat, paste0(opts$output,"/cnv.txt"),row.names = T,col.names = F, sep=",", quote=FALSE)


