#!/usr/bin/env Rscript
## Score SVCFit's hemizygous cellular-fraction estimates against the
## chrX simulation ground truth.
##
## WHY THIS DOES NOT CALL run_svcfit(). run_svcfit() reaches calc_svcf() through
## load_data() -> characterize_sv() -> annotate_cnv(), and that route needs
## heterozygous germline SNPs and FACETS segments. The simulation deliberately
## has NEITHER: a male chrX carries no het SNPs (that is the entire reason the
## hemizygous estimator exists), and FACETS is SNP-based, so 03_chrx_depth_
## segmentation.R replaces it with depth-derived cn_bar. Fabricating an
## anno_sv_cnv frame to satisfy the diploid path would mean inventing the very
## quantities under test.
##
## So this driver supplies VAF (from SVtyper) and cn_bar (from the segmentation)
## directly to SVCFit's EXPORTED hemizygous helpers. The arithmetic being
## validated is therefore the package's, not a reimplementation. What is
## reimplemented here is only the dispatch -- which helper applies to which row --
## transcribed from calc_svcf()'s hemizygous branch and stated explicitly below so
## it can be checked against that source.
##
## Usage: 13_chrx_score_svcfit.R [condition ...]     (default: all 45)

suppressMessages({
  library(SVCFit); library(dplyr); library(tidyr); library(stringr)
})

## WHICH SVCFit GOT LOADED. The scores depend on the SVCFit build, so the loaded
## package must export the hemizygous helpers, and its location is printed. Point
## R_LIBS at an SVCFit bedf5ef library (see tools/run_svcfit_self_evaluation.sh).
##
HEMI_API <- c("hemizygous_del_svcf", "hemizygous_dup_svcf", "resolve_hemizygous_svcf")
.missing <- setdiff(HEMI_API, ls("package:SVCFit"))
if (length(.missing)) {
  stop("SVCFit at ", find.package("SVCFit"), " exports no ",
       paste(.missing, collapse = ", "),
       ".\n  Point R_LIBS at an SVCFit bedf5ef library.", call. = FALSE)
}
cat("SVCFit:", find.package("SVCFit"), "\n")

CHRX_DIR   <- Sys.getenv("CHRX_DIR", unset = ".")

## Coverage tag. The truth table and the score file are one file per run, so both
## are coverage-named -- see the note in chrx_common.sh. REQUIRED, no default.
COV        <- Sys.getenv("COV", unset = "")
if (!nzchar(COV)) stop("COV is not set, and there is no default. Run e.g. COV=50 ...",
                       call. = FALSE)
ctag       <- paste0("c", COV)

## Replicate tag, matching chrx_common.sh: REP selects an independent draw of
## SHORtS read noise and each replicate has its own tree, so a replicate must not
## read another's calls or write over another's scores. REP=0 leaves every path
## exactly as it was. The truth table is NOT rep-tagged -- the planted SVs and
## their cellular fractions are identical across replicates, only read noise
## differs.
REP        <- Sys.getenv("REP", unset = "0")
rep_tag    <- if (REP == "0") "" else paste0("_rep", REP)

TRUTH      <- file.path(CHRX_DIR, "truth", paste0("chrx_sv_ground_truth_", ctag, ".tsv"))
CN_BAR_DIR <- file.path(CHRX_DIR, paste0("cn_bar10k_m", rep_tag))  # canonical join rule
CALLS_DIR  <- file.path(CHRX_DIR, paste0("calls",       rep_tag))
## CHRX_SCORE_ROOT redirects only the output, so a rescore with another SVCFit build reads the
## same calls and segmentation but never overwrites the accepted scoring_rep* tables.
SCORE_ROOT <- Sys.getenv("CHRX_SCORE_ROOT", unset = CHRX_DIR)
OUT_DIR    <- file.path(SCORE_ROOT, paste0("scoring",   rep_tag))

## Breakpoint match tolerance. Manta emits CIPOS=-100,100 on these calls, and the
## VISOR-planted coordinate and Manta's reported POS differ by one base even on
## exact calls, so an exact join would match almost nothing.
POS_TOL <- 250L

## A locus counts as copy-neutral when cn_bar is within this of 1 copy. Same
## tolerance hemizygous_del_svcf() uses for its depth cross-check.
CN_TOL <- 0.15

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- inputs ---------------------------------------------------------------

read_truth <- function() {
  t <- read.delim(TRUTH, stringsAsFactors = FALSE)
  stopifnot(all(c("experiment","condition","sv_id","CHROM","start","end",
                  "sv_class","cellular_fraction") %in% names(t)))
  t
}

## SVtyper VCF -> one row per called SV with the allele counts.
## RO/AO rather than AB: AB is rounded to 2 decimals in the FORMAT field, and the
## hemizygous helpers take counts, not a ratio.
read_calls <- function(exp_name, cond) {
  f <- file.path(CALLS_DIR, paste0(exp_name, "_", cond), "svtyp",
                 paste0("svt_", exp_name, "_", cond, ".vcf"))
  if (!file.exists(f)) return(NULL)
  ln <- readLines(f, warn = FALSE)
  ln <- ln[!startsWith(ln, "#")]
  if (!length(ln)) return(NULL)
  f9 <- do.call(rbind, strsplit(ln, "\t", fixed = TRUE))

  info <- f9[, 8]
  fmt  <- strsplit(f9[, 9], ":", fixed = TRUE)
  smp  <- strsplit(f9[, 10], ":", fixed = TRUE)
  getf <- function(i, key) {
    j <- match(key, fmt[[i]]); if (is.na(j) || j > length(smp[[i]])) NA_character_ else smp[[i]][j]
  }
  geti <- function(key) vapply(seq_along(fmt), function(i) suppressWarnings(as.numeric(getf(i, key))), numeric(1))

  data.frame(
    id       = f9[, 3],
    CHROM    = f9[, 1],
    POS      = as.integer(f9[, 2]),
    END      = as.integer(str_match(info, "(?:^|;)END=([0-9]+)")[, 2]),
    svtype   = str_match(info, "(?:^|;)SVTYPE=([A-Z]+)")[, 2],
    imprecise= grepl("(^|;)IMPRECISE(;|$)", info),
    RO       = geti("RO"),
    AO       = geti("AO"),
    stringsAsFactors = FALSE
  )
}

read_cn_bar <- function(exp_name, cond) {
  f <- file.path(CN_BAR_DIR, paste0(exp_name, "_", cond, "_cn_bar.csv"))
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE)
}

## ---- matching -------------------------------------------------------------

## Manta emits several records per event (the mate breakends of an inversion, and
## a PRECISE plus an IMPRECISE record for the same locus). Collapse to one row per
## truth SV, preferring a PRECISE record and then the one with the most evidence,
## so a single event cannot be scored twice.
match_calls_to_truth <- function(tr, cl) {
  if (is.null(cl) || !nrow(cl)) return(tr %>% mutate(RO = NA_real_, AO = NA_real_,
                                                     id = NA_character_, detected = FALSE))
  out <- lapply(seq_len(nrow(tr)), function(i) {
    hit <- which(cl$CHROM == tr$CHROM[i] &
                 abs(cl$POS - tr$start[i]) <= POS_TOL &
                 abs(cl$END - tr$end[i])   <= POS_TOL)
    if (!length(hit)) return(data.frame(RO = NA_real_, AO = NA_real_,
                                        id = NA_character_, detected = FALSE))
    h <- cl[hit, , drop = FALSE]
    h <- h[order(h$imprecise, -(h$AO + h$RO)), , drop = FALSE]
    data.frame(RO = h$RO[1], AO = h$AO[1], id = h$id[1], detected = TRUE)
  })
  cbind(tr, do.call(rbind, out))
}

## cn_bar is keyed on the SV's own start. Nearest within tolerance rather than an
## exact join, for the same one-base reason as above.
attach_cn_bar <- function(df, cb) {
  if (is.null(cb) || !nrow(cb)) return(df %>% mutate(cn_bar = NA_real_))
  df$cn_bar <- vapply(seq_len(nrow(df)), function(i) {
    d <- abs(cb$POS - df$start[i])
    if (min(d) > POS_TOL) NA_real_ else cb$cn_bar[which.min(d)]
  }, numeric(1))
  df
}

## BACKGROUND copy number: the copy state of the region the SV sits IN, measured
## from segments that do not overlap the SV itself.
##
## This is what decides which hemizygous form applies, and it cannot be taken from
## the SV's own cn_bar. A deletion drives its own span below 1 copy whether or not
## the surrounding region is amplified, so keying on the SV's cn_bar sends every
## deletion down the same branch. e4
## deletions sit inside amplifications (background ~1.4) and need
## resolve_hemizygous_svcf(); e1 and e2 deletions sit on unamplified background and
## need hemizygous_del_svcf().
##
## Taken as the length-weighted mean cn_bar of segments overlapping a flank of
## FLANK_BP on each side, excluding any part inside the SV span.
FLANK_BP <- 500000L

attach_background_cn <- function(df, seg) {
  if (is.null(seg) || !nrow(seg)) return(df %>% mutate(bg_cn = NA_real_))
  df$bg_cn <- vapply(seq_len(nrow(df)), function(i) {
    s <- df$start[i]; e <- df$end[i]
    lo <- s - FLANK_BP; hi <- e + FLANK_BP
    ov <- seg[seg$chrom == df$CHROM[i] & seg$end > lo & seg$start < hi, , drop = FALSE]
    if (!nrow(ov)) return(NA_real_)
    ## overlap with the flanks only, i.e. minus the SV's own span
    w <- pmax(0, pmin(ov$end, hi) - pmax(ov$start, lo)) -
         pmax(0, pmin(ov$end, e)  - pmax(ov$start, s))
    if (sum(w) <= 0) return(NA_real_)
    sum(ov$cn_bar * w) / sum(w)
  }, numeric(1))
  df
}

read_seg <- function(exp_name, cond) {
  f <- file.path(CHRX_DIR, "seg10k", paste0("chrx_seg_", exp_name, "_", cond, ".csv"))
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE)
}

## ---- dispatch, transcribed from calc_svcf()'s hemizygous branch -----------
##
##   tandem duplication  -> hemizygous_dup_svcf(cn_bar, r = 2); upper bound, since
##                          r is not identifiable from one locus
##   deletion            -> hemizygous_del_svcf(AO, RO, cn_bar), which picks depth
##                          when cn_bar < 1 and h2 otherwise. Deletions BYPASS the
##                          sign rule -- see the rationale in hemizygous.R.
##   copy-neutral other  -> SVCF = VAF                   (calc_svcf.R:92)
##   copy-altered other  -> resolve_hemizygous_svcf(AO, RO, cn_bar), which picks
##                          the ordering by the sign of the SV-first form
##
## bg_cn is still computed and written out, but no longer routes anything: the
## feasibility of the depth estimate turned out to be the better discriminator,
## and unlike background amplification it separates e2 from e4. Kept because it is
## the natural check on that claim.
score_rows <- function(d) {
  n <- nrow(d)
  svcf <- rep(NA_real_, n); status <- rep(NA_character_, n)
  bound <- rep(FALSE, n);   branch <- rep(NA_character_, n)

  vaf <- d$AO / (d$AO + d$RO)
  neutral <- is.finite(d$cn_bar) & abs(d$cn_bar - 1) <= CN_TOL
  amplified_bg <- is.finite(d$bg_cn) & d$bg_cn > 1 + CN_TOL

  is_dup <- d$sv_class == "tandem duplication"
  is_del <- d$sv_class == "deletion"

  for (i in seq_len(n)) {
    if (!isTRUE(d$detected[i]) || !is.finite(d$AO[i]) || !is.finite(d$RO[i])) {
      status[i] <- "not_detected"; branch[i] <- "none"; next
    }
    if (is_dup[i]) {
      r <- hemizygous_dup_svcf(cn_bar = d$cn_bar[i], r = 2)
      svcf[i] <- r$svcf; status[i] <- r$status; bound[i] <- isTRUE(r$is_upper_bound)
      branch[i] <- "dup_bound"
    } else if (is_del[i]) {
      ## bg_cn is the flanking copy number, already computed by attach_background_cn(). It supplies
      ## kappa for the CNV-first form; without it the helper falls back to its previous behaviour.
      r <- hemizygous_del_svcf(bpc = d$AO[i], bec = d$RO[i], cn_bar = d$cn_bar[i],
                               bg_cn = d$bg_cn[i])
      svcf[i] <- r$svcf; status[i] <- r$status
      branch[i] <- paste0("del_", if (!is.na(r$svcf_source)) r$svcf_source else "NA")
    } else if (isTRUE(neutral[i]) && !is_del[i]) {
      svcf[i] <- vaf[i]; status[i] <- if (is.finite(vaf[i])) "ok" else "zero_ref_depth"
      branch[i] <- "copy_neutral_vaf"
    } else if (is.finite(d$cn_bar[i])) {
      r <- resolve_hemizygous_svcf(bpc = d$AO[i], bec = d$RO[i], cn_bar = d$cn_bar[i])
      svcf[i] <- r$svcf; status[i] <- if (!is.null(r$status)) r$status else "ok"
      branch[i] <- paste0("resolved_", if (!is.null(r$ordering)) r$ordering else "NA")
    } else {
      status[i] <- "no_cn_bar"; branch[i] <- "none"
    }
  }
  d$vaf <- vaf; d$svcf <- svcf; d$svcf_status <- status
  d$is_upper_bound <- bound; d$branch <- branch
  d
}

## ---- run ------------------------------------------------------------------

truth <- read_truth()
conds <- unique(truth[, c("experiment", "condition")])
sel <- commandArgs(trailingOnly = TRUE)
if (length(sel)) conds <- conds[conds$condition %in% sel | paste0(conds$experiment,"_",conds$condition) %in% sel, ]

res <- lapply(seq_len(nrow(conds)), function(k) {
  e <- conds$experiment[k]; c_ <- conds$condition[k]
  tr <- truth[truth$experiment == e & truth$condition == c_, ]
  d <- match_calls_to_truth(tr, read_calls(e, c_))
  d <- attach_cn_bar(d, read_cn_bar(e, c_))
  d <- attach_background_cn(d, read_seg(e, c_))
  score_rows(d)
})
all <- bind_rows(res)
all$err <- all$svcf - all$cellular_fraction

out <- file.path(OUT_DIR, paste0("chrx_svcfit_scores_", ctag, ".tsv"))
write.table(all, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("wrote", out, "--", nrow(all), "truth SVs\n\n")

## Tandem-duplication upper bounds are INCLUDED in the error statistics, so this
## summary matches the bootstraps in 19/23. They are r=2 upper bounds, not point estimates -- the
## is_upper_bound column and the bounded count below keep that split visible. Only non-finite errors
## (undetected and infeasible-DUP rows) are dropped.
pt <- all %>% filter(is.finite(err))
cat("detection: ", sum(all$detected), "/", nrow(all),
    sprintf(" (%.1f%%)\n", 100*mean(all$detected)), sep = "")
cat("scored (point estimates): ", nrow(pt), "\n", sep = "")
cat(sprintf("mean |err| %.4f   median %.4f   p95 %.4f   within 0.05 %.1f%%\n\n",
            mean(abs(pt$err)), median(abs(pt$err)),
            quantile(abs(pt$err), .95), 100*mean(abs(pt$err) <= .05)))

cat("by experiment:\n")
print(pt %>% group_by(experiment) %>%
        summarise(n = n(), mean_abs_err = mean(abs(err)),
                  within_05 = 100*mean(abs(err) <= .05), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\nby sv_class:\n")
print(pt %>% group_by(sv_class) %>%
        summarise(n = n(), mean_abs_err = mean(abs(err)),
                  within_05 = 100*mean(abs(err) <= .05), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\nby branch:\n")
print(all %>% group_by(branch) %>% summarise(n = n(), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\nsvcf_status:\n")
print(as.data.frame(table(all$svcf_status, useNA = "ifany")), row.names = FALSE)

nb <- sum(all$is_upper_bound & is.finite(all$svcf))
cat("\nbounded (tandem duplication, r not identifiable): ", nb, " rows included above as upper bounds\n", sep="")
