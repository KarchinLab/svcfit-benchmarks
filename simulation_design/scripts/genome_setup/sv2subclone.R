# assign SV to clones
library(dplyr)
library(readr)
library(optparse)

option_list <- list(
  make_option(c("-a", "--hack1"),  
              type    = "character",  
              help    = "Path to hack1 file"),
  
  make_option(c("-b", "--hack2"),
              type    = "character",
              help    = "Path to hack2 file"),

 # make_option(c("-c", "--hack3"),
 #             type    = "character",
 #             help    = "Path to hack3 file"),

  make_option(c("-o", "--out"),
              type    = "character",
              help    = "output dir")
)

parser <- OptionParser(option_list = option_list,add_help_option = FALSE)
opts   <- parse_args(parser)

H1 <- read.delim(opts$hack1, header=FALSE)
H1_1 <- read.delim(opts$hack2, header=FALSE)
#aH1_1 <- read.delim(opts$hack3, header=FALSE)
colnames(H1)=c('CHROM','start','end','class','info','flank')
colnames(H1_1)=c('CHROM','start','end','class','info','flank')
#colnames(aH1_1)=c('CHROM','start','end','class','info','flank')

## h1 subclone
rows=round(nrow(H1)/3)
c1=H1[1:rows,]
tc2=H1[(rows+1):(2*rows),]
tc3=H1[(2*rows+1):nrow(H1),]
hom_c1=c1 %>% mutate(grp=row_number()%%2) %>% filter(grp==0) %>% select(-grp)
thom_c2=tc2 %>% mutate(grp=row_number()%%2) %>% filter(grp==0) %>% select(-grp)
thom_c3=tc3 %>% mutate(grp=row_number()%%2) %>% filter(grp==0) %>% select(-grp)

c2=rbind(tc2, c1) %>% arrange(CHROM,start)
c3=rbind(tc3, c1) %>% arrange(CHROM,start)
hom_c2=rbind(thom_c2, hom_c1) %>% arrange(CHROM,start)
hom_c3=rbind(thom_c3, hom_c1) %>% arrange(CHROM,start)


#rm(tc2, tc3, thom_c2, thom_c3, rows, dup1, dup2, dup3)


# h11 subclone
rows=round(nrow(H1_1)/3)
c11=H1_1[1:rows,]
tc22=H1_1[(rows+1):(2*rows),]
tc33=H1_1[(2*rows+1):nrow(H1_1),]

hom_c11=c11 %>% mutate(grp=row_number()%%2) %>% filter(grp==0) %>% select(-grp)
thom_c22=tc22 %>% mutate(grp=row_number()%%2) %>% filter(grp==0) %>% select(-grp)
thom_c33=tc33 %>% mutate(grp=row_number()%%2) %>% filter(grp==0) %>% select(-grp)

c22=rbind(tc22, c11) %>% arrange(CHROM,start)
c33=rbind(tc33, c11) %>% arrange(CHROM,start)
hom_c22=rbind(thom_c22, hom_c11) %>% arrange(CHROM,start)
hom_c33=rbind(thom_c33, hom_c11) %>% arrange(CHROM,start)

dup11=data.frame(CHROM='chr1', start=(min(c11$start)-1000), end=(max(c11$end)+1000), 
                class='tandem duplication', info=3, flank=0)
dup22=data.frame(CHROM='chr1', start=(min(tc22$start)-1000), end=(max(tc22$end)+1000), 
                class='tandem duplication', info=5, flank=0)
dup33=data.frame(CHROM='chr1', start=(min(tc33$start)-1000), end=(max(tc33$end)+1000), 
                class='tandem duplication', info=2, flank=0)
dup_c22=rbind(dup11,dup22)
dup_c33=rbind(dup11,dup33)


del11=data.frame(CHROM='chr1', start=(min(c11$start)-1000), end=(max(c11$end)+1000),
                class='deletion', info='None', flank=0)
del22=data.frame(CHROM='chr1', start=(min(tc22$start)-1000), end=(max(tc22$end)+1000),
                class='deletion', info='None', flank=0)
del33=data.frame(CHROM='chr1', start=(min(tc33$start)-1000), end=(max(tc33$end)+1000),
                class='deletion', info='None', flank=0)
del_c22=rbind(del11,del22)
del_c33=rbind(del11,del33)
#rm(tc22, tc33, thom_c22, thom_c33, rows)

# ## for sv after in-cis cnv
len=(dup11$end-dup11$start)*(dup11$info-1)
taft_c22=tc22 %>%
  mutate(start=start+len,
         end=end+len)
taft_c33=tc33 %>%
  mutate(start=start+len,
         end=end+len)
aft_c22=rbind(c11,taft_c22)
aft_c33=rbind(c11,taft_c33)
  

write_delim(c1, paste0(opts$out,"/c1.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(c2, paste0(opts$out,"/c2.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(c3, paste0(opts$out,"/c3.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(c11, paste0(opts$out,"/c11.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(c22, paste0(opts$out,"/c22.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(c33, paste0(opts$out,"/c33.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(hom_c1, paste0(opts$out,"/hom_c1.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(hom_c2, paste0(opts$out,"/hom_c2.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(hom_c3, paste0(opts$out,"/hom_c3.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(hom_c11, paste0(opts$out,"/hom_c11.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(hom_c22, paste0(opts$out,"/hom_c22.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(hom_c33, paste0(opts$out,"/hom_c33.bed"), quote = 'none', delim = "\t", col_names = F)
#write_delim(aft_c11, paste0(opts$out,"/aft_c11.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(aft_c22, paste0(opts$out,"/aft_c22.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(aft_c33, paste0(opts$out,"/aft_c33.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(dup_c22, paste0(opts$out,"/dup_c22.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(dup_c33, paste0(opts$out,"/dup_c33.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(del_c22, paste0(opts$out,"/del_c22.bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(del_c33, paste0(opts$out,"/del_c33.bed"), quote = 'none', delim = "\t", col_names = F)
