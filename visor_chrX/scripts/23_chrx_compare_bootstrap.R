#!/usr/bin/env Rscript
##
## Block bootstrap of the SVCFit-vs-SVclone comparison over the chrX replicates.
##
## Same resampling unit and machinery as 19_chrx_bootstrap.R: replicates are drawn with replacement,
## B = 10000, percentile CIs, matching the submitted autosomal figure
## (the autosomal benchmark). Read 19's header before quoting anything from here; every
## caveat there applies unchanged, in particular:
##
##   THE REPLICATES SHARE ONE MATCHED NORMAL. 02's normal mode uses a hardcoded seed that REP never
##   reaches, so a per-replicate normal would be byte-identical. cn_bar is a tumour-to-normal depth
##   ratio, so the normal's read noise is a constant rather than a resampled quantity, and these
##   intervals resample TUMOUR read noise only. Intervals on depth-dependent branches are narrower
##   than a full resampling would give. That is a property of the simulation, and it applies to both
##   tools equally, so it does not bias the comparison -- but it must travel with the numbers.
##
## THE PAIRED DIFFERENCE IS THE STATISTIC THAT MATTERS. Both tools see the same replicate: the same
## BAM, the same calls, the same segments. Resampling replicates and taking the difference WITHIN
## each draw keeps that pairing, so replicate-to-replicate noise common to both tools cancels. Two
## independent per-tool intervals would not cancel it, and overlapping marginal intervals would be
## read as "no difference" when the paired difference is unambiguous. Both are reported; the paired
## one is the answer.
##
## HEMIZYGOUS TANDEM DUPLICATIONS ARE INCLUDED, as upper bounds. SVCFit cannot
## point-identify a tandem duplication's cellular fraction from one locus -- r is not identifiable --
## so hemizygous_dup_svcf() returns the r=2 upper bound. We COUNT the upper bounds as SVCFit's
## estimate in every stratum, including ALL, so the comparison covers the same SVs for both tools.
## The only remaining filter is `is.finite(err)`, which still drops the infeasible DUP rows (cn_bar
## below one copy -> svcf NA). READ THE DUP ROW WITH THIS IN MIND: SVCFit's value there is an upper
## bound on the cellular fraction, not a point estimate, and depth-based bounds at low CF are poor,
## so a worse SVCFit DUP error is expected and is not evidence the estimator is miscalibrated.
##
## Usage: COV=50 CHRX_DIR=... 23_chrx_compare_bootstrap.R [--out FILE]

suppressMessages({ library(dplyr); library(tidyr) })

CHRX_DIR <- Sys.getenv("CHRX_DIR", unset = ".")
COV      <- Sys.getenv("COV", unset = "")
if (!nzchar(COV)) stop("COV is not set, and there is no default. Run e.g. COV=50 ...", call. = FALSE)
B    <- as.integer(Sys.getenv("B", unset = "10000"))
SEED <- as.integer(Sys.getenv("BOOT_SEED", unset = "1"))
ctag <- paste0("c", COV)

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 2 && args[1] == "--out") args[2] else
        file.path(CHRX_DIR, "scoring", paste0("chrx_compare_bootstrap_", ctag, ".tsv"))

## ---- load both tools, per replicate -------------------------------------------------------------
dirs <- sort(Sys.glob(file.path(CHRX_DIR, "scoring_rep*")))
load_tool <- function(file, tool) {
  out <- list()
  for (d in dirs) {
    f <- file.path(d, file)
    if (!file.exists(f)) next
    rep <- as.integer(sub(".*scoring_rep([0-9]+).*", "\\1", d))
    x <- read.delim(f, stringsAsFactors = FALSE)
    x$replicate <- rep; x$tool <- tool
    out[[length(out) + 1]] <- x
  }
  bind_rows(out)
}
a <- load_tool(sprintf("chrx_svcfit_scores_%s.tsv",  ctag), "SVCFit")
b <- load_tool(sprintf("chrx_svclone_scores_%s.tsv", ctag), "SVclone")
if (!nrow(a)) stop("no SVCFit replicate tables found under ", CHRX_DIR, call. = FALSE)
if (!nrow(b)) stop("no SVclone replicate tables found under ", CHRX_DIR, call. = FALSE)

reps <- sort(intersect(unique(a$replicate), unique(b$replicate)))
cat(sprintf("replicates with BOTH tools: %d  (%s)\n", length(reps), paste(reps, collapse = ", ")))
only_a <- setdiff(unique(a$replicate), reps); only_b <- setdiff(unique(b$replicate), reps)
## A replicate scored for one tool and not the other cannot enter a paired comparison, and silently
## keeping it for the marginal intervals would make the two columns rest on different data.
if (length(only_a)) cat(sprintf("  DROPPED, SVCFit only:  %s\n", paste(only_a, collapse = ", ")))
if (length(only_b)) cat(sprintf("  DROPPED, SVclone only: %s\n", paste(only_b, collapse = ", ")))
if (length(reps) < 10)
  warning("only ", length(reps), " paired replicates: a percentile interval from this few is a ",
          "description of spread, not a confidence interval", call. = FALSE)
a <- a[a$replicate %in% reps, ]; b <- b[b$replicate %in% reps, ]

## ---- per-replicate, per-stratum statistics ------------------------------------------------------
## Each tool filtered by its OWN rule; see the header. SVCFit's tandem-duplication upper bounds are
## now kept (only non-finite errors -- undetected and infeasible DUP rows -- are dropped).
pa <- a %>% filter(is.finite(err))
pb <- b %>% filter(detected, is.finite(err))

stats_of <- function(pt, all) {
  data.frame(detection    = 100 * mean(all$detected),
             within_05    = if (nrow(pt)) 100 * mean(abs(pt$err) <= 0.05) else NA_real_,
             mean_abs_err = if (nrow(pt)) mean(abs(pt$err)) else NA_real_,
             n            = nrow(pt))
}
per_rep <- function(pt, all, tool) {
  bind_rows(
    all %>% group_by(replicate) %>% group_modify(~stats_of(pt[pt$replicate == .y$replicate, ], .x)) %>%
      mutate(stratum = "ALL", .before = 1),
    all %>% group_by(replicate, experiment) %>%
      group_modify(~stats_of(pt[pt$replicate == .y$replicate & pt$experiment == .y$experiment, ], .x)) %>%
      rename(stratum = experiment) %>% relocate(stratum),
    all %>% group_by(replicate, sv_class) %>%
      group_modify(~stats_of(pt[pt$replicate == .y$replicate & pt$sv_class == .y$sv_class, ], .x)) %>%
      rename(stratum = sv_class) %>% relocate(stratum)
  ) %>% mutate(tool = tool)
}
ra <- per_rep(pa, a, "SVCFit"); rb <- per_rep(pb, b, "SVclone")

## ---- bootstrap ----------------------------------------------------------------------------------
set.seed(SEED)
strata  <- sort(unique(c(ra$stratum, rb$stratum)))
metrics <- c("detection", "within_05", "mean_abs_err")
draws   <- replicate(B, sample(reps, length(reps), replace = TRUE), simplify = FALSE)

pick <- function(d, s, m) { v <- d[d$stratum == s, ]; setNames(v[[m]], v$replicate) }
res <- list()
for (s in strata) for (m in metrics) {
  va <- pick(ra, s, m); vb <- pick(rb, s, m)
  if (!length(va) || !length(vb)) next
  ## Draw once per bootstrap replicate and index BOTH tools with it: that is what keeps the pairing.
  bs <- vapply(draws, function(ix) {
    ka <- va[as.character(ix)]; kb <- vb[as.character(ix)]
    c(mean(ka, na.rm = TRUE), mean(kb, na.rm = TRUE),
      mean(ka, na.rm = TRUE) - mean(kb, na.rm = TRUE))
  }, numeric(3))
  q <- function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  res[[length(res) + 1]] <- data.frame(
    stratum = s, metric = m, n_rep = length(reps),
    svcfit  = mean(va, na.rm = TRUE),  svcfit_lo  = q(bs[1, ])[1], svcfit_hi  = q(bs[1, ])[2],
    svclone = mean(vb, na.rm = TRUE),  svclone_lo = q(bs[2, ])[1], svclone_hi = q(bs[2, ])[2],
    diff    = mean(va, na.rm = TRUE) - mean(vb, na.rm = TRUE),
    diff_lo = q(bs[3, ])[1], diff_hi = q(bs[3, ])[2])
}
out <- bind_rows(res)

## ---- report -------------------------------------------------------------------------------------
cat(sprintf("\nB = %d draws, seed %d, resampling %d replicates with replacement\n",
            B, SEED, length(reps)))
cat("diff = SVCFit - SVclone, paired within each draw. A CI excluding 0 is a real difference.\n")
for (s in strata) {
  d <- out[out$stratum == s, ]
  if (!nrow(d)) next
  cat(sprintf("\n%-20s\n", s))
  for (m in metrics) {
    r <- d[d$metric == m, ]
    if (!nrow(r)) next
    fmt <- if (m == "mean_abs_err") "%.4f" else "%.1f%%"
    star <- if (is.finite(r$diff_lo) && is.finite(r$diff_hi) && (r$diff_lo > 0 | r$diff_hi < 0)) " *" else ""
    cat(sprintf("  %-13s SVCFit %s   SVclone %s   diff %s [%s, %s]%s\n", m,
                sprintf(fmt, r$svcfit), sprintf(fmt, r$svclone), sprintf(fmt, r$diff),
                sprintf(fmt, r$diff_lo), sprintf(fmt, r$diff_hi), star))
  }
}
cat("\n* = 95% CI on the paired difference excludes zero\n")
cat("\nNOTE: for the tandem duplication stratum SVCFit's value is an r=2 UPPER BOUND on the\n")
cat("cellular fraction, not a point estimate (r is not identifiable from one locus); SVclone\n")
cat("reports a point estimate. The bounds are included in every stratum, ALL among them, so both\n")
cat("tools cover the same SVs -- but a worse SVCFit DUP error reflects the bound, not miscalibration.\n")

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.table(out, OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nwrote %s  (%d rows)\n", OUT, nrow(out)))
