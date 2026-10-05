library(dplyr)
library(readr)
library(optparse)

option_list <- list(
  make_option(c("-s", "--sample"),
              type = "character",
              help = "sample_name"),
  make_option(c("-o", "--o_path"),
              type = "character",
              help = "where to save output"),
  make_option(c("-v", "--sv_file"),
              type = "character",
              help = "where is SV file")
)

parser <- OptionParser(option_list = option_list,add_help_option = FALSE)
opts   <- parse_args(parser)

sv <- read.table(opts$sv_file, quote="\"")
colnames(sv) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','normal','tumor')

positions <- sv %>%
  mutate(chr2=ifelse(grepl('BND',ID), gsub(".*(chr\\d+):.*",'\\1',ALT), CHROM),
         pos2=ifelse(grepl('BND',ID), gsub(".*:(\\d+).*",'\\1',ALT) ,gsub(".*END=(\\d+);.*",'\\1', INFO)))

left=positions %>%
  select(CHROM, POS)%>%
  mutate(lpos=as.integer(POS)-500,
         rpos=as.integer(POS)+500)%>%
  select(-POS)

right=positions %>%
  select(chr2, pos2)%>%
  mutate(lpos=as.integer(pos2)-500,
         rpos=as.integer(pos2)+500,
         CHROM=chr2)%>%
  select(-chr2, -pos2)

all=rbind(left,right)%>%
  group_by(CHROM, lpos, rpos)%>%
  distinct(, .keep_all = T)

write_delim(all, paste0(opts$o_path, "/",opts$sample, ".bed"), quote = 'none', delim = "\t", col_names = F)
