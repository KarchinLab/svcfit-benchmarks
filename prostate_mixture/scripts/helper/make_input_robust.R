#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(dplyr); library(optparse) })
opts <- parse_args(OptionParser(option_list=list(
  make_option(c("-f","--file"),type="character"), make_option(c("-o","--output"),type="character"),
  make_option(c("-c","--cnv"),type="character"), make_option(c("-s","--smp"),type="character"))))
required <- c("file","output","cnv","smp")
missing <- required[vapply(required,function(x)is.null(opts[[x]])||!nzchar(opts[[x]]),logical(1))]
if(length(missing))stop("Missing options: ",paste(missing,collapse=", "))
dir.create(opts$output,recursive=TRUE,showWarnings=FALSE)

sv <- read.delim(opts$file,header=FALSE,comment.char="#",quote="",fill=TRUE,
                 stringsAsFactors=FALSE,check.names=FALSE)
if(ncol(sv)!=11L)stop("Expected exactly 11 VCF columns; found ",ncol(sv))
names(sv)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","normal","tumor")
info <- function(x,key){m<-regexec(paste0("(?:^|;)",key,"=([^;]+)"),x,perl=TRUE);z<-regmatches(x,m);vapply(z,function(y)if(length(y)>=2)y[2]else NA_character_,character(1))}
class<-info(sv$INFO,"SVTYPE"); is_bnd<-class=="BND"|grepl("[\\[\\]]",sv$ALT,perl=TRUE)
m<-regexec("[\\[\\]]([^:\\[\\]]+):(\\d+)[\\[\\]]",sv$ALT,perl=TRUE);z<-regmatches(sv$ALT,m)
mate_chr<-vapply(z,function(x)if(length(x)>=3)x[2]else NA_character_,character(1));mate_pos<-vapply(z,function(x)if(length(x)>=3)x[3]else NA_character_,character(1))
pos1<-suppressWarnings(as.integer(sv$POS));chr2<-ifelse(is_bnd,mate_chr,sv$CHROM)
pos2<-suppressWarnings(as.integer(ifelse(is_bnd,mate_pos,info(sv$INFO,"END"))))
keep<-sv$FILTER=="PASS"&class!="INS";bad<-keep&(is.na(class)|is.na(pos1)|is.na(chr2)|is.na(pos2))
if(any(bad))stop("Malformed retained coordinates: ",paste(sv$ID[bad],collapse=", "))
simple<-data.frame(chr1=sv$CHROM[keep],pos1=pos1[keep],chr2=chr2[keep],pos2=pos2[keep])
inv<-class[keep]=="INV";if(any(inv)){x<-simple$pos1[inv];simple$pos1[inv]<-simple$pos2[inv];simple$pos2[inv]<-x}
write.table(simple,file.path(opts$output,paste0("svc_",opts$smp,"_simple.txt")),row.names=FALSE,sep="\t",quote=FALSE)

e<-new.env(parent=emptyenv());load(file.path(opts$cnv,paste0(opts$smp,".RData")),envir=e)
if(!exists("fit",envir=e,inherits=FALSE))stop("FACETS fit object missing")
fit<-get("fit",envir=e);pur<-as.numeric(fit$purity);ploidy<-as.numeric(fit$ploidy)
if(length(pur)!=1L||!is.finite(pur))stop("FACETS purity unavailable; truth fallback forbidden")
if(length(ploidy)!=1L||!is.finite(ploidy))stop("FACETS ploidy unavailable")
write.table(data.frame(sample=opts$smp,purity=pur,ploidy=ploidy),file.path(opts$output,paste0("spp_",opts$smp,".txt")),row.names=FALSE,sep="\t",quote=FALSE)
facet<-read.delim(file.path(opts$cnv,paste0(opts$smp,".bed")))
need<-c("chrom","start","end","tcn.em","lcn.em");if(length(setdiff(need,names(facet))))stop("Invalid FACETS bed")
cnv<-facet%>%mutate(chrom=as.character(chrom),chrom=ifelse(chrom=="23","X",ifelse(chrom=="24","Y",chrom)),n_major=2,n_minor=1,t_total=tcn.em,t_minor=ifelse(is.na(lcn.em),0L,lcn.em))%>%select(chrom,start,end,n_major,n_minor,t_total,t_minor)
write.table(cnv,file.path(opts$output,"cnv.txt"),row.names=TRUE,col.names=FALSE,sep=",",quote=FALSE)
message("Wrote ",nrow(simple)," retained SVs; FACETS purity=",pur)
