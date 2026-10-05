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

  make_option(c("-o", "--out"),
              type    = "character",
              help    = "Path to output file"),
  
  make_option(c("-n", "--nsv"),
              type    = "integer",
              help    = "num of sv needed"),
  
  make_option(c("-z", "--size"),
              type    = "integer",
              help    = "sv size")
)

parser <- OptionParser(option_list = option_list,add_help_option = FALSE)
opts   <- parse_args(parser)

print(paste0("The SV bed is ",opts$hack1))
print(paste0("The bed for transloc is ",opts$hack2))
print(paste0("output file is located ",opts$out))

sample = ifelse(grepl('HACk11', opts$hack1), "H1_1", "H1")

c1 <- read.delim(opts$hack1, header=FALSE)[5:(opts$nsv+4),]
c2 <- read.delim(opts$hack2, header=FALSE)

colnames(c1) <- c('chr','start','end','sv','info','flank')
colnames(c2) <- c('chr','start4','end4','sv','info','flank')

ntrans=sum(c1$sv=="reciprocal translocation")

## create various trans between chr1 and chr2
trans <- c1 %>%
  filter(sv=="reciprocal translocation")%>%
  cbind(c2[5:(5+ntrans-1),2:3])%>%
  mutate(id=row_number()%%3,
         sv=case_when(
           id==1 ~ "reciprocal translocation",
           id==2 ~ "translocation copy-paste",
           id==0 ~ "translocation cut-paste"),
         info=gsub("chr\\d","chr2", info),
         info=gsub("h2","h1", info),
         info=ifelse(id==1, info,gsub(":\\w+$","", info)))%>%
  rowwise()%>%
  mutate(info=gsub(":(\\d+):",paste0(":",start4,":"),info))%>%
  ungroup()%>%
  select(-start4,-end4,-id)

## create dup with various copies
if (sum(c1$sv=="tandem duplication")>0){
	dup <- c1 %>%
		filter(sv=="tandem duplication")%>%
		mutate(info = sample(2:5,sum(c1$sv=="tandem duplication"), replace=T))
} else {
	dup=data.frame()
}

## replace the old dup and trans with new trans on 2 chrom and dup with various copies 
## make SV length vary as well
hack <- c1%>%
  filter(sv!="reciprocal translocation")%>%
  filter(sv!="tandem duplication")%>%
  rbind(dup)%>%
  rbind(trans)%>%
  arrange(chr,start)%>%
  mutate(end=start+round(rnorm(nrow(.), mean = opts$size, sd = 10000)))

## make every other sv homozygous
hom=hack %>%
  mutate(type=row_number()%%2)%>%
  filter(type==0)%>%
  select(-type)

write_delim(hack, paste0(opts$out, "/",sample,".bed"), quote = 'none', delim = "\t", col_names = F)
write_delim(hom, paste0(opts$out, "/hom_",sample,".bed"), quote = 'none', delim = "\t", col_names = F)
