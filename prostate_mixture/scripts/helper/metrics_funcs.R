library(ggplot2)
library(gridExtra)
library(data.table)
library(plyr)
library(GenomicRanges)
library(gtools)
library(RColorBrewer)
library(stringr)


attach_001_ground_truth_svs <- function(wd_truth, ccf_all, ff=15) {
    bm_svs <- read.table(paste(wd_truth, "/bM.vcf", sep = ''), sep='\t', stringsAsFactors = F, header = F)
    gm_svs <- read.table(paste(wd_truth, "/gM.vcf", sep = ''), sep='\t', stringsAsFactors = F, header = F)
    colnames(bm_svs) <- colnames(gm_svs) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','normal','tumor')
    
    bm_svs <- bm_svs %>% 
      filter(FILTER=='PASS') %>% 
      mutate(pp=ifelse(grepl('BND', ID), 
                         str_extract(ALT, "(?<=:)\\d+"), 
                         gsub(".*END=(\\d+);.*","\\1", INFO)),
             chr2=ifelse(grepl('BND', ID), 
                         str_extract(ALT, "(?<=\\[|\\])[^:]+"), 
                         CHROM),
             POS=as.numeric(POS),
             pos2=as.numeric(pp))%>%
      ## chrX is included, matching parse_truth_sv. chr2 also gets 'X' so chrX BNDs mating to X and
      ## the X-ends of BNDs survive; Y and decoy mates are excluded because they are not in the list.
      filter(CHROM %in% c(1:22, 'X'), chr2 %in% c(1:22, 'X'))
    
    gm_svs <- gm_svs %>% 
      filter(FILTER=='PASS') %>% 
      mutate(pp=ifelse(grepl('BND', ID), 
                         str_extract(ALT, "(?<=:)\\d+"), 
                         gsub(".*END=(\\d+);.*","\\1", INFO)),
             chr2=ifelse(grepl('BND', ID), 
                         str_extract(ALT, "(?<=\\[|\\])[^:]+"), 
                         CHROM),
             POS=as.numeric(POS),
             pos2=as.numeric(pp))%>%
      ## chrX is included, matching parse_truth_sv. chr2 also gets 'X' so chrX BNDs mating to X and
      ## the X-ends of BNDs survive; Y and decoy mates are excluded because they are not in the list.
      filter(CHROM %in% c(1:22, 'X'), chr2 %in% c(1:22, 'X'))
    
    # check SVs that don't match exactly
    bm_svs$id <- paste(bm_svs$CHROM, bm_svs$POS, 
                       bm_svs$chr2, bm_svs$pos2, sep=':')
    gm_svs$id <- paste(gm_svs$CHROM, gm_svs$POS,
                       gm_svs$chr2, gm_svs$pos2, sep=':')
    
    both <- c()
    for (i in 1:nrow(bm_svs)) {
      bm_chr1 <- bm_svs[i,'CHROM']
      bm_chr2 <- bm_svs[i,'chr2']
      bm_pos1 <- bm_svs[i,'POS']
      bm_pos2 <- bm_svs[i,'pos2']
      for (j in 1:nrow(gm_svs)) {
        gm_chr2 <- gm_svs[j,'chr2']
        gm_chr1 <- gm_svs[j,'CHROM']
        gm_pos1 <- gm_svs[j,'POS']
        gm_pos2 <- gm_svs[j,'pos2']
        if (gm_chr1 == bm_chr1 & gm_chr2 == bm_chr2) {
          if (abs(gm_pos1 - bm_pos1) < ff & abs(gm_pos2 - bm_pos2) < ff) {
            both <- c(both, bm_svs$id[i], gm_svs$id[j])
          }
        }
        
      }
    }
    both <- unique(both)
    bm_uniq <- bm_svs$id[(!bm_svs$id%in%both)]
    gm_uniq <- gm_svs$id[(!gm_svs$id%in%both)]
    
    ccf_all$sample <- NA
    ccf_all[ccf_all$sv%in%bm_uniq, 'sample'] <- 'bM_only'
    ccf_all[ccf_all$sv%in%gm_uniq, 'sample'] <- 'gM_only'
    ccf_all[ccf_all$sv%in%both, 'sample'] <- 'shared'
    ccf_all <- ccf_all[!is.na(ccf_all$sample), ]
    
    return(ccf_all)
}

attach_001_ground_truth_snvs <- function(snv_wd, ccfs) {
    bm_snvs <- suppressWarnings(fread(paste(snv_wd, "001bM_SS_mut_merged.call.stats.txt", sep = ''), 
                                      sep='\t', stringsAsFactors = F, header = T))
    gm_snvs <- suppressWarnings(fread(paste(snv_wd, "001gM_SS_mut_merged.call.stats.txt", sep = ''), 
                                      sep='\t', stringsAsFactors = F, header = T))
    bm <- bm_snvs$Row.names
    gm <- gm_snvs$Row.names
    
    both <- intersect(bm, gm)
    gm_uniq <- gm[!gm%in%bm]
    bm_uniq <- bm[!bm%in%gm]
    
    ccfs[ccfs$sv%in%bm_uniq, 'sample'] <- 'bM_only'
    ccfs[ccfs$sv%in%gm_uniq, 'sample'] <- 'gM_only'
    ccfs[ccfs$sv%in%both, 'sample'] <- 'shared'
    ccfs[is.na(ccfs$sample), 'sample'] <- 'unmatched'
    unmatched <- ccfs[ccfs$sample == 'unmatched', ]
    ccfs <- ccfs[ccfs$sample!='unmatched', ]
    
    return(ccfs)
}

# Metrics for 4- and 5-cluster solutions; calc_metrics_3clus below handles 3 clusters.
calc_metrics_hiclus <- function(x, mix, nclus=c(4,5), type='sv', truth='best') {
    method <- NA
    if(type == 'svs') {
        mean_mult_error <- mean(c(x$true_cn1-x$ccube_mult1, x$true_cn2-x$ccube_mult2))
        mean_ccf_error <- mean(c(x$best_ccf1-x$ccube_ccf1, x$best_ccf2-x$ccube_ccf2))
        ccfs <- as.numeric(apply(x[,c('ccube_ccf1', 'ccube_ccf2')], 1, mean))
        best_ccf <- as.numeric(apply(x[,c('best_ccf1', 'best_ccf2')], 1, mean))
        if(truth=='ground'){best_ccf <- x$true_ccf}
        method <- 'ccube_SV'
        cutoff <- cutoff_svs
        clus_ccf <- x$ccube_ccf_mean
    } else if(type == 'snvs') {
        mean_mult_error <- mean(x$true_cn - x$ccube_mult)
        mean_ccf_error <- mean(x$best_ccf - x$ccube_ccf)
        ccfs <- x$ccube_ccf
        best_ccf <- x$best_ccf
        if(truth=='ground'){best_ccf <- x$true_ccf}
        method <- 'ccube_SNV'
        cutoff <- cutoff_snvs
        clus_ccf <- x$ccube_ccf_mean
    } else {
        mean_mult_error <- mean(x$true_cn - x$multiplicity)
        mean_ccf_error <- mean(x$best_ccf - x$ccf)
        ccfs <- x$ccf
        best_ccf <- x$best_ccf
        if(truth=='ground'){best_ccf <- x$true_ccf}
        method <- 'pyclone'
        cutoff <- cutoff_snvs
        clus_ccf <- x$cluster_ccf
    }
    if(nclus==4){truth <- c(0.2, 0.4, 0.6, 1)}else{truth <- c(0.2, 0.4, 0.6, 0.8, 1)}
    clus_num_error <- nclus - length(unique(clus_ccf))
    clus_ccf_error <- get_clus_ccf_error(clus_ccf, truth)
    is_subclonal_truth <- best_ccf < cutoff
    is_subclonal_ccube <- ccfs < cutoff
    is_subclonal_sensitivity <- sum(is_subclonal_truth & is_subclonal_ccube) / sum(is_subclonal_truth)
    is_subclonal_specificity <- 1 - sum(!is_subclonal_truth & is_subclonal_ccube) / nrow(x)
    met <- data.frame(mix, mean_mult_error, mean_ccf_error, clus_num_error,
                      is_subclonal_sensitivity, is_subclonal_specificity,
                      clus_ccf_error, method=method)
    return(met)
}

calc_metrics_3clus <- function(x, mix, method, type='sv') {
    tmp <- as.numeric(strsplit(mix, '-')[[1]])
    tmp <- tmp[order(tmp, decreasing = T)]
    major <- tmp[1]; minor <- tmp[2]
    
    ccfs <- NULL
    if(type=='pyc') {
        mean_mult_error <- mean(c(x$true_cn - x$multiplicity))
        mean_ccf_error <- mean(c(x$true_ccf - x$ccf))
        clus_num_error <- 3 - length(unique(x$cluster_ccf))
        clus_ccf_error <- get_clus_ccf_error(x$cluster_ccf, c(1, major, minor))
        is_subclonal_truth <- x$true_ccf < 1
        is_subclonal <- x$ccf < cutoff_snvs
    } else if(type=='sv') {
        x$mean_ccf <- apply(data.frame(ccf1=x$ccube_ccf1, ccf2=x$ccube_ccf2), 1, mean)
        mean_mult_error <- mean(c(x$true_cn1-x$ccube_mult1, x$true_cn2-x$ccube_mult1))
        mean_ccf_error <- mean(c(x$true_ccf-x$ccube_ccf1, x$true_ccf-x$ccube_ccf2))
        clus_num_error <- 3 - length(unique(x$ccube_ccf_mean))
        clus_ccf_error <- get_clus_ccf_error(x$ccube_ccf_mean, c(1, major, minor))
        is_subclonal_truth <- x$true_ccf < 1
        is_subclonal <- apply(x[,c('ccube_ccf1', 'ccube_ccf2')], 1, mean) < cutoff_svs
    } else if (type=='snv') {
        mean_mult_error <- mean(x$true_cn - x$ccube_mult)
        mean_ccf_error <- mean(x$true_ccf - x$ccube_ccf)
        clus_num_error <- 3 - length(unique(x$ccube_ccf_mean))
        clus_ccf_error <- get_clus_ccf_error(x$ccube_ccf_mean, c(1, major, minor))
        is_subclonal_truth <- x$true_ccf < 1
        is_subclonal <- x$ccube_ccf< cutoff_snvs
    }
    is_subclonal_sensitivity <- sum(is_subclonal_truth & is_subclonal) / sum(is_subclonal_truth)
    is_subclonal_specificity <- 1 - sum(!is_subclonal_truth & is_subclonal) / nrow(x)
    # precision <- sum(is_subclonal_truth & is_subclonal) / length(is_subclonal_truth)
    # f1_score <- 2 * ((precision * is_subclonal_sensitivity)/(precision + is_subclonal_sensitivity))
    return(data.frame(mix, mean_mult_error, mean_ccf_error, clus_num_error,
                      is_subclonal_sensitivity, is_subclonal_specificity,
                      clus_ccf_error, method=method))
    
}

get_clus_ccf_error <- function(cms, truth) {
    # lowest ccf cluster = match with minor fraction
    # highest ccf clsuter = match with clonal
    # average the ccfs of remaining clusters, match with major
    cm <- as.numeric(names(table(cms)))
    cm <- sort(cm, decreasing=T)
    truth <- sort(truth, decreasing=T)
    clus_ccf_error <- NULL
    for(i in 1:length(truth)) {
        if(length(cm)==0){break}
        clus_ccf_error <- c(clus_ccf_error, truth[1] - cm[1])
        cm <- cm[!cm%in%cm[1]]
        cm <- rev(cm)
        truth <- rev(truth)
    }
    return(mean(clus_ccf_error))
}