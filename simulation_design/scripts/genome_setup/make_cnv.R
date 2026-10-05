library(dplyr)
library(readr)
library(optparse)

option_list <- list(
  make_option(c("-n", "--nsv"),
              type    = "integer",
              help    = "num of sv needed"),
  make_option(c("-p", "--in_path"),
              type    = "character",
              help    = "input file path"),
  make_option(c("-o", "--out_path"),
              type    = "character",
              help    = "input file path")
)

parser <- OptionParser(option_list = option_list,add_help_option = FALSE)
opts   <- parse_args(parser)

half_num=opts$nsv/2

## this uses modified SV with new trans and dup
hack = read.delim(paste0(opts$in_path,"/H1_1.bed"), header=FALSE)
colnames(hack) <- c('chr','start','end','sv','info','flank')

## create overlapping dup
dup_over <- hack %>%
  mutate(type=rep(1:2, each=half_num),
         svlen=ifelse(sv=="deletion", end-start, 0))%>%
  group_by(type)%>%
  mutate(tlen=sum(svlen),
         start=min(start)-2000,
         end=max(end)+2000,
         sv="tandem duplication",
         info=ifelse(type==1, "5", "2"),
         flank=0)%>%
  ungroup()%>%
  distinct(type, .keep_all = T)%>%
  mutate(end=ifelse(type==1, end-tlen, end-sum(tlen)),
         start=ifelse(type==2, start-tlen[1], start))%>%
  select(-type, -svlen, -tlen)

## create overlapping del
del_over = dup_over %>%
  mutate(sv="deletion",
         info="None")
del1=del_over[1,]
del2=del_over[2,]
del1$end=del2$end

### now adjusting sv pposition for cnv_sv (dup)
## foe del, never same phase with SV, no need to correct position
dup1bp=(dup_over$end[1]-dup_over$start[1])*(as.integer(dup_over$info[1])-1)
hack[(half_num+1):nrow(hack),]$start = hack[(half_num+1):nrow(hack),]$start + dup1bp
hack[(half_num+1):nrow(hack),]$end = hack[(half_num+1):nrow(hack),]$end + dup1bp
hack[1:half_num,]$start = hack[1:half_num,]$start + dup1bp
hack[1:half_num,]$end = hack[1:half_num,]$end + dup1bp

hom=hack %>%
  mutate(type=row_number()%%2)%>%
  filter(type==0)%>%
  select(-type)

write_delim(dup_over, paste0(opts$out_path,"/dup_over.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(del1,paste0(opts$out_path, "/del1.bed"), quote = 'none', delim = "\t", col_names = F)
#write_delim(del2, "del2.bed", quote = 'none', delim = "\t", col_names = F)
write_delim(hack, paste0(opts$out_path,"/aH1_1.bed"), quote = 'none', delim = "\t", col_names = F)
