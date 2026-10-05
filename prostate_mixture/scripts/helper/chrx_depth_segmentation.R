#!/usr/bin/env Rscript
# chrX copy-number segmentation from read depth, for male subjects.
#
# HARD PREREQUISITE for the hemizygous path (MATH.md changes 3 and 4). It emits ONE quantity:
#
#   cn_bar  the mean copy number per cell of a chrX segment
#
# It deliberately does NOT emit c or f_CNV. Depth cannot separate them: cn_bar = 1 + f(c-1) is one
# equation in two unknowns, and a male X has no heterozygous SNPs to supply the allelic-imbalance
# equation FACETS uses on the autosomes, and nothing downstream needs c or f_CNV.
#
# Neither is otherwise available. SVCFit computes cellular fractions for structural variants, not
# for copy-number segments. FACETS fits a male X with no heterozygous SNPs and returns one
# whole-chromosome segment with an inflated tcn.em, which this pipeline rejects. And the containing
# structural variant's own SVCFit estimate is circular, because it is hemizygous and copy-altered too.
#
# Why depth works here: the matched normal is also single-copy on chrX, so after normalising each
# BAM by its own autosomal depth the tumour-to-normal ratio at a chrX segment is the mean copy number
# per cell. With cells carrying the change at c copies and the rest at 1,
#
#     cn_bar = f*c + (1 - f)*1 = 1 + f*(c - 1)   =>   f = (cn_bar - 1) / (c - 1)
#
#   Rscript chrx_depth_segmentation.R <tumor.bam> <normal.bam> <psi_sample> <out.csv> [bin_kb]
#
# psi_sample is the mean autosomal copy number per cell IN THE SAMPLE, i.e. including normal
# contamination: psi_sample = purity * psi_tumour + (1 - purity) * 2, with psi_tumour the
# length-weighted mean of the FACETS autosomal tcn.em. FACETS is trustworthy here: the autosomes have
# heterozygous SNPs, so this is the fit it is designed for, unlike its chrX fit which this pipeline
# rejects. Across this cohort psi_sample/2 runs from 0.87 to 1.72, so it is not optional.
#
# Requires samtools on PATH and the DNAcopy package.

suppressPackageStartupMessages(library(DNAcopy))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: chrx_depth_segmentation.R <tumor.bam> <normal.bam> <psi_sample> <out.csv> [bin_kb]")
tbam <- args[1]; nbam <- args[2]; psi_sample <- as.numeric(args[3]); out <- args[4]
bin_kb <- if (length(args) >= 5) as.numeric(args[5]) else 100
if (!is.finite(psi_sample) || psi_sample <= 0) stop("psi_sample must be a positive number")

## ASSEMBLY-DEPENDENT, and overridable by environment variable. The prostate
## mixture benchmark is GRCh37/hs37d5: its X contig is named "X", not "chrX", and is
## 155,270,560 bp rather than hg38's 156,040,895. Left as hardcoded hg38 values these would have
## failed in the worst possible way -- `samtools depth -r chrX:...` on a GRCh37 BAM returns no rows,
## every bin reads zero depth, and the script emits a segmentation of an empty chromosome rather
## than an error.
##
## Defaults are the previous hardcoded values EXACTLY, so the COMBAT cohort path is unchanged and
## no existing caller passes these. Set all three together; mixing an hg38 length with GRCh37 names
## silently truncates or overruns the contig.
##   CHRX_CONTIG    contig name as it appears in the BAM           (default chrX)
##   CHRX_LEN       its length in bp                               (default 156040895, hg38)
##   CHRX_AUTO_REF  comma-separated autosomal 2-copy baseline regions, same naming as the contig
CHRX_CONTIG <- Sys.getenv("CHRX_CONTIG", unset = "chrX")
CHRX_LEN    <- as.numeric(Sys.getenv("CHRX_LEN", unset = "156040895"))
AUTO_REF    <- strsplit(Sys.getenv("CHRX_AUTO_REF",
                 unset = "chr1:50000000-52000000,chr2:100000000-102000000,chr7:20000000-22000000,chr12:60000000-62000000"),
                 ",")[[1]]
if (!is.finite(CHRX_LEN) || CHRX_LEN <= 0) stop("CHRX_LEN must be a positive number", call. = FALSE)
## A region naming mismatch is the failure this guards: it produces zero depth everywhere rather
## than an error, so it is checked here instead of being discovered in the output.
if (any(grepl("^chr", AUTO_REF)) != grepl("^chr", CHRX_CONTIG))
  stop("CHRX_CONTIG ('", CHRX_CONTIG, "') and CHRX_AUTO_REF disagree on the 'chr' prefix; ",
       "one of them does not match the BAM", call. = FALSE)

mean_depth <- function(bam, region, mapq = 20) {
  cmd <- sprintf("samtools depth -a -Q %d -r %s %s | awk '{s+=$3; n++} END{if(n) print s/n; else print 0}'",
                 mapq, region, shQuote(bam))
  as.numeric(system(cmd, intern = TRUE))
}

## autosomal baseline: what 2 copies looks like in each BAM
auto <- function(bam) median(vapply(AUTO_REF, function(r) mean_depth(bam, r), numeric(1)))
t_auto <- auto(tbam); n_auto <- auto(nbam)
cat(sprintf("autosomal mean depth: tumour %.1f, normal %.1f\n", t_auto, n_auto))
if (t_auto <= 0 || n_auto <= 0) stop("zero autosomal depth; check the BAM paths")

## binned chrX depth
starts <- seq(1, CHRX_LEN - bin_kb * 1000, by = bin_kb * 1000)
cat(sprintf("binning %s at %d kb, %d bins\n", CHRX_CONTIG, bin_kb, length(starts)))
td <- nd <- numeric(length(starts))
for (i in seq_along(starts)) {
  r <- sprintf("%s:%d-%d", CHRX_CONTIG, starts[i], starts[i] + bin_kb * 1000 - 1)
  td[i] <- mean_depth(tbam, r); nd[i] <- mean_depth(nbam, r)
}
## Zero depth in EVERY bin of both BAMs means the contig name does not match, not that the
## chromosome is absent -- samtools returns nothing for an unknown region and says nothing about it.
if (all(td == 0) && all(nd == 0))
  stop("every bin of ", CHRX_CONTIG, " read zero depth in both BAMs. The contig name almost ",
       "certainly does not match the BAM header; check CHRX_CONTIG.", call. = FALSE)
## Zero in the TUMOUR alone is a different fault and needs its own message. The contig name is fine
## -- the normal found reads on it -- so this BAM genuinely lacks the chromosome. Without this the
## run continues to cn_bar = 0, log2(0) = -Inf, and dies inside DNAcopy with "only 0's may be mixed
## with negative subscripts", which names neither the chromosome nor the BAM.
if (all(td == 0))
  stop("every bin of ", CHRX_CONTIG, " read zero depth in the TUMOUR but not the normal, so the ",
       "contig name is right and this BAM has no ", CHRX_CONTIG, ": ", tbam, call. = FALSE)

## normalise each by its own autosomal (2-copy) level, then take the ratio.
## normal chrX is 1 copy, so the ratio is the tumour's mean copy number per cell.
##
## * psi_sample/2 RECOVERS ABSOLUTE COPY NUMBER. The tumour's autosomes sit at
## psi_sample copies, not 2, so
##
##     td / t_auto            = cn_bar / psi_sample
##     nd / n_auto            = 1 / 2
##  => (td/t_auto)/(nd/n_auto) = 2 * cn_bar / psi_sample
##
## and the factor of psi_sample/2 turns that back into cn_bar.
keep <- nd > 0 & td > 0
cn_bar <- (td[keep] / t_auto) / (nd[keep] / n_auto) * (psi_sample / 2)
pos <- starts[keep]
cat(sprintf("%d informative bins, median cn_bar %.2f\n", sum(keep), median(cn_bar)))

## segment
## SEGMENTATION PARAMETERS, overridable and defaulting to DNAcopy's own.
##
## Defaults are unchanged deliberately: on the simulation every setting tried
## gives an IDENTICAL answer -- 1 segment on a flat chrX, 7 on one with real
## amplifications, with those amplifications (2.58, 1.32, 1.70) preserved in
## every case. The knobs do nothing on clean data, so changing the default on
## simulation evidence would be changing it on no evidence at all.
##
## They exist because real chrX is not clean data: 83694, a sample with flat
## depth and one SV, segments into 41 pieces spanning 0.83-0.96. That spread is
## too small to be copy number and too systematic to be bin noise -- GC and
## mappability differences between a tumour and normal library will do it.
##
## Measure before setting these: bin a sample
## once and segments it every way, which is the cheap way round: binning is two
## samtools depth calls per bin, segmentation is instant.
##
##   SEG_ALPHA      significance for a split          (DNAcopy default 0.01)
##   SEG_UNDO_SD    merge segments closer than this   (unset = no merging)
##                  many SDs apart; 1.5-2.0 is the usual range
##   SEG_MIN_WIDTH  minimum bins per segment          (DNAcopy default 2)
seg_alpha     <- as.numeric(Sys.getenv("SEG_ALPHA", "0.01"))
seg_undo_sd   <- Sys.getenv("SEG_UNDO_SD", "")
seg_min_width <- as.integer(Sys.getenv("SEG_MIN_WIDTH", "2"))

cna <- CNA(log2(cn_bar), rep(CHRX_CONTIG, length(pos)), pos,
           data.type = "logratio", sampleid = CHRX_CONTIG)
seg_args <- list(smooth.CNA(cna), verbose = 0,
                 alpha = seg_alpha, min.width = seg_min_width)
if (nzchar(seg_undo_sd)) {
  seg_args$undo.splits <- "sdundo"
  seg_args$undo.SD     <- as.numeric(seg_undo_sd)
}
cat(sprintf("segmenting: alpha=%g min.width=%d undo.SD=%s\n",
            seg_alpha, seg_min_width, if (nzchar(seg_undo_sd)) seg_undo_sd else "none"))
seg <- do.call(segment, seg_args)$output

## per-segment mean copy number. No c, no f: see the header.
## chrom is CHRX_CONTIG, not the literal "chrX": the segment-to-SV join matches on this name, so
## labelling a GRCh37 run's output "chrX" while its SVs sit on "X" would match nothing and report
## every SV as having no depth rather than failing.
res <- data.frame(chrom = CHRX_CONTIG,
                  start = seg$loc.start, end = seg$loc.end,
                  n_bins = seg$num.mark,
                  cn_bar = round(2^seg$seg.mean, 4),
                  stringsAsFactors = FALSE)

write.csv(res, out, row.names = FALSE)
print(res, row.names = FALSE)
cat(sprintf("\nwritten to %s\n", out))
cat("cn_bar is the mean copy number per cell, already scaled by psi_sample, and is the final\n")
cat("product of this step. The corrected hemizygous forms take it directly, so it is never rounded\n")
cat("to an integer copy number and there is no c/f solve downstream (see CORRECTION-c-vs-cnbar).\n")
cat("Sanity check: any chrX segment that is genuinely at germline copy number should read near 1.0.\n")
cat("If none does, either chrX is wholly altered in this sample or psi_sample is wrong.\n")
