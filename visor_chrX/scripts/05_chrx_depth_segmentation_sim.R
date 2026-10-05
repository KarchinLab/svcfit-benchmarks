#!/usr/bin/env Rscript
# chrX copy-number segmentation for the SIMULATED hemizygous BAMs.
#
# Derived from prostate_mixture/scripts/helper/chrx_depth_segmentation.R with one change for the
# simulation, marked [SIM] below.
#
# [SIM] AUTO_REF. The general script baselines on chr1/chr2/chr7/chr12. The simulated BAMs carry
#   ONLY chr22 and chrX (01_chrx_hack.sh builds h1 = chr22+chrX, h2 = chr22), so the general
#   AUTO_REF returns zero depth and the script stops at "zero autosomal depth". The four windows
#   used here are resources/beds/depth_windows_auto.bed, already verified 0.00% N.
#
#   THIS IS NOT COSMETIC: `samtools depth -a` scores assembly N at zero, but
#   SHORtS simulates reads only for non-N bases. chr22 is 22.9% N (acrocentric p-arm), so an
#   N-containing baseline window reads low and inflates every cn_bar on the chromosome. The
#   binned chrX ratio is immune — the same positions are N in tumour and normal, so the N
#   fraction cancels in td/nd — but the AUTO_REF baseline is a bare mean and does not cancel.
#
# psi_sample is applied as in the general script: cn_bar = ratio * psi_sample/2. For this
#   simulation chr22 is copy-neutral diploid in every clone by construction, so psi_sample = 2
#   and the factor is exactly 1.
#
#   Usage: Rscript 05_chrx_depth_segmentation_sim.R <tumor.bam> <normal.bam> <psi_sample> <out.csv> [bin_kb]
#
# Requires samtools on PATH and the DNAcopy package (both present in conda env `visor`).

suppressPackageStartupMessages(library(DNAcopy))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: 05_chrx_depth_segmentation_sim.R <tumor.bam> <normal.bam> <psi_sample> <out.csv> [bin_kb]")
tbam <- args[1]; nbam <- args[2]; psi_sample <- as.numeric(args[3]); out <- args[4]
bin_kb <- if (length(args) >= 5) as.numeric(args[5]) else 100
if (!is.finite(psi_sample) || psi_sample <= 0) stop("psi_sample must be a positive number")

CHRX_LEN <- 156040895                      # hg38
# [SIM 1] chr22 baseline, from resources/beds/depth_windows_auto.bed (verified N-free)
AUTO_REF <- c("chr22:20000000-22000000", "chr22:28000000-30000000",
              "chr22:36000000-38000000", "chr22:46000000-48000000")

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
cat(sprintf("binning chrX at %d kb, %d bins\n", bin_kb, length(starts)))
td <- nd <- numeric(length(starts))
for (i in seq_along(starts)) {
  r <- sprintf("chrX:%d-%d", starts[i], starts[i] + bin_kb * 1000 - 1)
  td[i] <- mean_depth(tbam, r); nd[i] <- mean_depth(nbam, r)
}

## normalise each by its own autosomal (2-copy) level, then take the ratio.
## normal chrX is 1 copy, so the ratio is the tumour's mean copy number per cell.
## [SIM 2] * psi_sample/2 recovers absolute mean copy number; exactly 1 here (psi_sample = 2).
keep <- nd > 0 & td > 0
cn_bar <- (td[keep] / t_auto) / (nd[keep] / n_auto) * (psi_sample / 2)
pos <- starts[keep]
cat(sprintf("%d informative bins, median cn_bar %.2f\n", sum(keep), median(cn_bar)))

## segment
cna <- CNA(log2(cn_bar), rep("chrX", length(pos)), pos, data.type = "logratio", sampleid = "chrX")
seg <- segment(smooth.CNA(cna), verbose = 0)$output

## per-segment mean copy number. No c, no f: see the header.
res <- data.frame(chrom = "chrX",
                  start = seg$loc.start, end = seg$loc.end,
                  n_bins = seg$num.mark,
                  cn_bar = round(2^seg$seg.mean, 4),
                  stringsAsFactors = FALSE)

write.csv(res, out, row.names = FALSE)
print(res, row.names = FALSE)
cat(sprintf("\nwritten to %s\n", out))
