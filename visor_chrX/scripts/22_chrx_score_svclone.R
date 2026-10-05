#!/usr/bin/env Rscript
##
## Score SVclone against the chrX ground truth, on the same terms as 13_chrx_score_svcfit.R.
##
## The point of this file is comparability. Wherever a choice exists it is made the way 13 makes it,
## so the two score tables can be put side by side without an argument about method:
##   - the same truth table, truth/chrx_sv_ground_truth_c${COV}.tsv
##   - the same breakpoint match tolerance, POS_TOL = 250 (Manta emits CIPOS=-100,100 and the
##     VISOR-planted coordinate differs from Manta's POS by a base even on exact calls)
##   - detection, within 0.05, and mean |err| computed identically
##
## SVCLONE'S OWN ESTIMATE IS USED, NOT A RECOMPUTED ONE. An SV's estimate is the proportion of the
## cluster SVclone assigned it to: most_likely_assignment in cluster_certainty.txt, looked up in
## subclonal_structure.txt. That is SVclone's answer, produced by SVclone's own read counting and
## clustering, and it is already on the truth table's scale -- a single clonal cluster at 10% purity
## reports proportion 0.1, which is the cellular fraction, not a cancer-cell fraction needing
## conversion.
##
## WHY NOT THE WRAPPER'S CCF PATH. The autosomal analysis
## recomputes CCF through get_true_sc_cn_by_side() -> MapVaf2CcfPyClone(). That works on autosomes
## and does not transfer to a hemizygous chromosome. SVclone's read-counting convention runs about
## twice the naive allelic fraction on BOTH (autosomal c50p10m10: vaf1 median 0.101 where a clonal
## het SV at 10% purity implies 0.05; chrX: 0.156 where it implies 0.10). On a diploid locus
## total_cn = 2, the multiplicity search may pick m = 2, and the factor divides out to give CCF 1.01.
## On a hemizygous locus total_cn = 1, so `ms <- 1:cn[1]` offers only m = 1, the factor has nowhere
## to go, and CCF comes out ~1.56 instead of 1. The compensation is structural to the diploid case
## and has no hemizygous analogue. Using SVclone's own output sidesteps it entirely, and also means
## this script needs neither ccube nor PROSTATE_R_DIR.
##
## Usage: 22_chrx_score_svclone.R [--scale-check] [condition ...]     (default: all)

suppressMessages(library(dplyr))

CHRX_DIR <- Sys.getenv("CHRX_DIR", unset = ".")
COV      <- Sys.getenv("COV", unset = "")
if (!nzchar(COV)) stop("COV is not set, and there is no default. Run e.g. COV=50 ...", call. = FALSE)
REP      <- as.integer(Sys.getenv("REP", unset = "0"))
rep_tag  <- if (REP > 0) paste0("_rep", REP) else ""

args        <- commandArgs(trailingOnly = TRUE)
scale_check <- "--scale-check" %in% args
want_conds  <- setdiff(args, "--scale-check")

TRUTH   <- file.path(CHRX_DIR, "truth", sprintf("chrx_sv_ground_truth_c%s.tsv", COV))
SVC_DIR <- file.path(CHRX_DIR, paste0("svclone", rep_tag))
OUT_DIR <- file.path(CHRX_DIR, paste0("scoring", rep_tag))
POS_TOL <- 250L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

truth <- read.delim(TRUTH, stringsAsFactors = FALSE)
stopifnot(all(c("experiment","condition","sv_id","CHROM","start","end",
                "sv_class","cellular_fraction") %in% names(truth)))

conds <- unique(paste(truth$experiment, truth$condition, sep = "_"))
if (length(want_conds)) conds <- intersect(conds, want_conds)
if (!length(conds)) stop("no conditions selected", call. = FALSE)

## ---- one condition -----------------------------------------------------------------------------
score_one <- function(name) {
  exp_name <- sub("_.*", "", name); cond <- sub("^[^_]*_", "", name)
  tr <- truth[truth$experiment == exp_name & truth$condition == cond, , drop = FALSE]
  if (!nrow(tr)) return(NULL)
  tr$svcf_svclone <- NA_real_; tr$detected <- FALSE
  tr$cluster <- NA_integer_; tr$tool <- "svclone"

  base <- file.path(SVC_DIR, name, name, "ccube_out")
  cert <- file.path(base, paste0(name, "_cluster_certainty.txt"))
  strc <- file.path(base, paste0(name, "_subclonal_structure.txt"))
  pp_p <- file.path(SVC_DIR, name, name, "purity_ploidy.txt")
  ## A condition SVclone could not fit is not dropped: every truth SV in it counts as undetected,
  ## exactly as it would if the caller had found nothing. Dropping it instead would quietly raise
  ## SVclone's accuracy by removing the cases it failed on.
  if (!all(file.exists(cert, strc))) {
    message(sprintf("[%s] no SVclone fit -- %d truth SV(s) counted as undetected", name, nrow(tr)))
    tr$purity <- NA_real_
    return(tr)
  }

  cc <- read.delim(cert, stringsAsFactors = FALSE)
  ss <- read.delim(strc, stringsAsFactors = FALSE)
  pur <- if (file.exists(pp_p)) as.numeric(read.delim(pp_p, stringsAsFactors = FALSE)$purity[1]) else NA_real_

  ## SVclone's estimate for an SV is the proportion of the cluster it was assigned to.
  cc$prop <- ss$proportion[match(cc$most_likely_assignment, ss$cluster)]

  cs <- sub("^chr", "", as.character(cc$chr1))
  ct <- sub("^chr", "", as.character(cc$chr2))
  p1 <- suppressWarnings(as.numeric(cc$pos1)); p2 <- suppressWarnings(as.numeric(cc$pos2))

  for (i in seq_len(nrow(tr))) {
    tc <- sub("^chr", "", tr$CHROM[i])
    hit <- which(cs == tc & ct == tc &
                 abs(p1 - tr$start[i]) <= POS_TOL & abs(p2 - tr$end[i]) <= POS_TOL)
    if (!length(hit))
      ## make_input.R swaps pos1/pos2 for inversions on the way in, so a reversed match is the same
      ## call, not a second one.
      hit <- which(cs == tc & ct == tc &
                   abs(p2 - tr$start[i]) <= POS_TOL & abs(p1 - tr$end[i]) <= POS_TOL)
    if (length(hit)) {
      j <- hit[1]
      tr$detected[i]     <- TRUE
      tr$svcf_svclone[i] <- cc$prop[j]
      tr$cluster[i]      <- cc$most_likely_assignment[j]
    }
  }
  tr$err    <- tr$svcf_svclone - tr$cellular_fraction
  tr$purity <- pur
  tr
}

res <- bind_rows(lapply(conds, function(cn) tryCatch(score_one(cn), error = function(e) {
  warning(sprintf("[%s] %s", cn, conditionMessage(e))); NULL })))
if (!nrow(res)) stop("no conditions scored", call. = FALSE)

## ---- scale check -------------------------------------------------------------------------------
## A clonal SV sits in every tumour cell, so its cellular_fraction equals the purity. If SVclone's
## reported proportion is on the same scale, the ratio is 1. This is checked rather than assumed,
## because a scale error would leave every number below looking plausible and being wrong.
if (scale_check) {
  cl <- res %>% filter(detected, is.finite(svcf_svclone), is.finite(purity),
                       abs(cellular_fraction - purity) < 1e-6)
  cat("\n--- scale check: clonal SVs, where estimate/purity should be 1 ---\n")
  if (!nrow(cl)) {
    cat("  no clonal SVs matched; scale unverified\n")
  } else {
    r <- cl$svcf_svclone / cl$purity
    cat(sprintf("  n = %d   median ratio = %.3f   (median estimate %.4f vs truth %.4f)\n",
                nrow(cl), median(r), median(cl$svcf_svclone), median(cl$cellular_fraction)))
    cat(if (abs(median(r) - 1) <= 0.15) "  scale OK\n" else
        "  ** RATIO IS NOT 1 -- the scale is wrong and the numbers below are meaningless **\n")
  }
}

## ---- report, in 13's shape ---------------------------------------------------------------------
pt <- res %>% filter(detected, is.finite(err))
cat(sprintf("\ndetection: %d/%d (%.1f%%)\n", sum(res$detected), nrow(res), 100 * mean(res$detected)))
cat(sprintf("scored (point estimates): %d\n", nrow(pt)))
if (nrow(pt)) {
  cat(sprintf("mean |err| %.4f   median %.4f   p95 %.4f   within 0.05 %.1f%%\n",
              mean(abs(pt$err)), median(abs(pt$err)),
              quantile(abs(pt$err), 0.95), 100 * mean(abs(pt$err) <= 0.05)))
  stat <- function(d) d %>% summarise(n = n(), mean_abs_err = mean(abs(err)),
                                      within_05 = 100 * mean(abs(err) <= 0.05))
  cat("\nby experiment:\n"); print(as.data.frame(pt %>% group_by(experiment) %>% stat()), row.names = FALSE)
  cat("\nby sv_class:\n");   print(as.data.frame(pt %>% group_by(sv_class)   %>% stat()), row.names = FALSE)
}

out <- file.path(OUT_DIR, sprintf("chrx_svclone_scores_c%s.tsv", COV))
write.table(res, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nRESULT: wrote %s (%d rows)\n", out, nrow(res)))
