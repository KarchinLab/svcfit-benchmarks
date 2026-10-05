#!/usr/bin/env Rscript

# Build SVclone inputs using FACETS purity and robust BND/END parsing. Ground
# truth purity is never used as an estimate; it is accepted only as a label.
suppressPackageStartupMessages({ library(dplyr); library(optparse) })

options <- list(
  make_option(c("-f", "--file"), type="character"),
  make_option(c("-p", "--truth-label"), type="character", dest="truth_label"),
  make_option(c("-o", "--output"), type="character"),
  make_option(c("-c", "--cnv"), type="character"),
  make_option(c("-s", "--smp"), type="character"),
  make_option(c("-r", "--rep"), type="character"),
  make_option(c("-e", "--exp"), type="character"),
  make_option(c("-t", "--purity-table"), type="character", dest="purity_table")
)
opts <- parse_args(OptionParser(option_list=options))
required <- c("file","output","cnv","smp","rep","exp","purity_table")
missing <- required[vapply(required, function(x) is.null(opts[[x]]) || !nzchar(opts[[x]]), logical(1))]
if (length(missing)) stop("Missing options: ", paste(missing, collapse=", "))
dir.create(opts$output, recursive=TRUE, showWarnings=FALSE)

status_path <- file.path(opts$output, paste0("facets_status_", opts$smp, ".txt"))
write_status <- function(status, purity=NA_real_, note="") {
  write.table(data.frame(sample=opts$smp, arm="fair", facets_status=status,
    purity_used=purity, truth_purity_label=opts$truth_label %||% "", note=note),
    status_path, row.names=FALSE, sep="\t", quote=FALSE)
}
`%||%` <- function(x,y) if (is.null(x)) y else x

fac <- read.csv(opts$purity_table, stringsAsFactors=FALSE)
row <- fac[fac$rep == opts$rep & fac$exp == opts$exp & fac$cond == opts$smp,,drop=FALSE]
if (nrow(row) != 1L) stop("Expected one purity row; found ", nrow(row))
purity <- suppressWarnings(as.numeric(row$purity[1])); ploidy <- suppressWarnings(as.numeric(row$ploidy[1]))
if (!is.finite(purity)) {
  write_status("failed_no_purity", note="FACETS purity unavailable; no truth substitution")
  message("FAIR_SKIP_NO_PURITY")
  quit(status=3)
}
if (!is.finite(ploidy)) stop("FACETS ploidy is unavailable")

sv <- read.delim(opts$file, header=FALSE, comment.char="#", quote="", fill=TRUE,
                 stringsAsFactors=FALSE, check.names=FALSE)
if (ncol(sv) != 11L) stop("Expected exactly 11 VCF columns; found ", ncol(sv))
names(sv) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","normal","tumor")
extract_info <- function(x,key) {
  m <- regexec(paste0("(?:^|;)",key,"=([^;]+)"),x,perl=TRUE); z <- regmatches(x,m)
  vapply(z,function(y) if(length(y)>=2L)y[2] else NA_character_,character(1))
}
class <- extract_info(sv$INFO,"SVTYPE")
is_bnd <- class == "BND" | grepl("[\\[\\]]",sv$ALT,perl=TRUE)
mm <- regexec("[\\[\\]]([^:\\[\\]]+):(\\d+)[\\[\\]]",sv$ALT,perl=TRUE); mp <- regmatches(sv$ALT,mm)
mate_chr <- vapply(mp,function(x) if(length(x)>=3L)x[2] else NA_character_,character(1))
mate_pos <- vapply(mp,function(x) if(length(x)>=3L)x[3] else NA_character_,character(1))
pos1 <- suppressWarnings(as.integer(sv$POS)); chr2 <- ifelse(is_bnd,mate_chr,sv$CHROM)
pos2 <- suppressWarnings(as.integer(ifelse(is_bnd,mate_pos,extract_info(sv$INFO,"END"))))
keep <- sv$FILTER == "PASS" & class != "INS"
bad <- keep & (is.na(class)|is.na(pos1)|is.na(chr2)|is.na(pos2))
if (any(bad)) stop("Malformed retained coordinates: ",paste(sv$ID[bad],collapse=", "))
simple <- data.frame(chr1=sv$CHROM[keep],pos1=pos1[keep],chr2=chr2[keep],pos2=pos2[keep])
inv <- class[keep] == "INV"; if(any(inv)){x<-simple$pos1[inv];simple$pos1[inv]<-simple$pos2[inv];simple$pos2[inv]<-x}
write.table(simple,file.path(opts$output,paste0("svc_",opts$smp,"_simple.txt")),row.names=FALSE,sep="\t",quote=FALSE)
write.table(data.frame(sample=opts$smp,purity=purity,ploidy=ploidy),
  file.path(opts$output,paste0("spp_",opts$smp,".txt")),row.names=FALSE,sep="\t",quote=FALSE)
facet <- read.delim(file.path(opts$cnv,paste0(opts$smp,".bed")))
need <- c("chrom","start","end","tcn.em","lcn.em"); if(length(setdiff(need,names(facet)))) stop("Invalid FACETS bed")
cnv <- facet %>% mutate(chrom=paste0("chr",chrom),n_major=2,n_minor=1,t_total=tcn.em,
  t_minor=ifelse(is.na(lcn.em),0L,lcn.em)) %>% select(chrom,start,end,n_major,n_minor,t_total,t_minor)
write.table(cnv,file.path(opts$output,"cnv.txt"),row.names=TRUE,col.names=FALSE,sep=",",quote=FALSE)
write_status("ok",purity,"FACETS purity; native SVclone multiplicity")
message("Wrote ",nrow(simple)," retained SVs; FACETS purity=",purity)
