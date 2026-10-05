#!/usr/bin/env Rscript
## Block bootstrap over the chrX simulation replicates.
##
##   COV=50 CHRX_DIR=<visor_chrX> Rscript 19_chrx_bootstrap.R [--out FILE]
##
## THE RESAMPLING UNIT IS THE REPLICATE, and that is the whole point. Each
## replicate is an independent draw of VISOR SHORtS read noise, which is what the
## submitted autosomal figure resamples (the autosomal benchmark, B =
## 10000, percentile CIs) and the reason replicate support was added at all.
##
## Resampling SVs instead would be wrong and would look fine. SVs within a
## condition share one BAM, one segmentation and one purity, so resampling them
## measures SV-to-SV heterogeneity, not read noise, and returns tighter intervals
## than the submitted method for a reason that has nothing to do with estimator
## quality. If anyone ever reports SV-level intervals they must be labelled as
## such and never placed beside the submitted figure.
##
## WHAT THESE INTERVALS DO NOT COVER. The replicates share ONE matched normal --
## 02's normal mode uses a hardcoded seed that REP never reaches, so a
## per-replicate normal would be byte-identical. cn_bar is a tumour-to-normal
## depth ratio, so the normal's read noise is a constant across every replicate
## and is NOT resampled here. Intervals on the cn_bar-dependent branches (e2, e4,
## depth-routed deletions) are therefore narrower than a full resampling would
## give. Copy-neutral inversions, where SVCF = VAF and no depth enters, are
## unaffected. Say so wherever these numbers appear.
##
## Metric definitions are taken from 13_chrx_score_svcfit.R and not re-invented.
## HEMIZYGOUS TANDEM DUPLICATIONS ARE INCLUDED as upper bounds: they count
## toward the error statistics in every stratum, ALL among them, as in 13 and
## 23. The only filter is `is.finite(err)`, which still drops undetected and infeasible-DUP rows.
## n_bounded is still reported so the bound-vs-point-estimate split stays visible; the DUP value is
## an r=2 upper bound, not a point estimate, so a larger DUP error is the bound, not miscalibration.

suppressMessages({ library(dplyr); library(tidyr) })

CHRX_DIR <- Sys.getenv("CHRX_DIR", unset = ".")
COV      <- Sys.getenv("COV", unset = "")
if (!nzchar(COV)) stop("COV is not set, and there is no default. Run e.g. COV=50 ...", call. = FALSE)
B        <- as.integer(Sys.getenv("B", unset = "10000"))
SEED     <- as.integer(Sys.getenv("BOOT_SEED", unset = "1"))
ctag     <- paste0("c", COV)

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 2 && args[1] == "--out") {
  args[2]
} else {
  file.path(CHRX_DIR, "scoring", paste0("chrx_bootstrap_", ctag, ".tsv"))
}

## ---- load every replicate ---------------------------------------------------

dirs <- sort(Sys.glob(file.path(CHRX_DIR, paste0("scoring_rep*"))))
files <- file.path(dirs, paste0("chrx_svcfit_scores_", ctag, ".tsv"))
files <- files[file.exists(files) & file.size(files) > 0]
if (!length(files)) stop("no replicate score tables found under ", CHRX_DIR, call. = FALSE)

reps <- sub(".*scoring_rep([0-9]+).*", "\\1", files)
cat(sprintf("replicates found: %d  (%s)\n", length(files), paste(sort(as.integer(reps)), collapse = ", ")))

## A percentile interval from a handful of replicates is not a confidence
## interval, it is a rounding of the smallest and largest value seen. Warn rather
## than refuse, so partial runs can still be inspected.
if (length(files) < 10)
  warning("only ", length(files), " replicates: percentile CIs from fewer than ~10 ",
          "are dominated by the extremes and should not be quoted.", call. = FALSE)

dat <- bind_rows(lapply(seq_along(files), function(i) {
  d <- read.delim(files[i], stringsAsFactors = FALSE)
  d$replicate <- as.integer(reps[i])
  d
}))

need <- c("experiment", "sv_class", "detected", "err", "is_upper_bound", "svcf")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("score tables lack: ", paste(miss, collapse = ", "), call. = FALSE)

dat$detected       <- as.logical(dat$detected)
dat$is_upper_bound <- as.logical(dat$is_upper_bound)

## ---- per-replicate statistics, matching 13 ---------------------------------

stats_of <- function(d) {
  ## Upper bounds count as estimates; only non-finite errors are dropped.
  pt <- d[is.finite(d$err), ]
  data.frame(
    n            = nrow(d),
    detection    = 100 * mean(d$detected),
    n_scored     = nrow(pt),
    mean_abs_err = if (nrow(pt)) mean(abs(pt$err)) else NA_real_,
    within_05    = if (nrow(pt)) 100 * mean(abs(pt$err) <= 0.05) else NA_real_,
    n_bounded    = sum(d$is_upper_bound & is.finite(d$svcf))
  )
}

per_rep <- bind_rows(
  dat %>% group_by(replicate) %>% group_modify(~stats_of(.x)) %>% mutate(stratum = "ALL", .before = 1),
  dat %>% group_by(replicate, experiment) %>% group_modify(~stats_of(.x)) %>%
    rename(stratum = experiment),
  dat %>% group_by(replicate, sv_class) %>% group_modify(~stats_of(.x)) %>%
    rename(stratum = sv_class)
) %>% ungroup()

## ---- the bootstrap ----------------------------------------------------------
## Resample REPLICATES with replacement, average their per-replicate statistics,
## repeat B times, take percentiles. Same shape as the submitted boot_f1b().

set.seed(SEED)
rep_ids <- sort(unique(per_rep$replicate))
n_rep   <- length(rep_ids)
metrics <- c("detection", "within_05", "mean_abs_err")

boot_one <- function(strat_df) {
  draws <- replicate(B, {
    s <- sample(rep_ids, n_rep, replace = TRUE)
    idx <- match(s, strat_df$replicate)
    colMeans(strat_df[idx, metrics, drop = FALSE], na.rm = TRUE)
  })
  as.data.frame(t(apply(draws, 1, function(v)
    c(ci_lo = unname(quantile(v, 0.025, na.rm = TRUE)),
      ci_hi = unname(quantile(v, 0.975, na.rm = TRUE))))))
}

res <- bind_rows(lapply(unique(per_rep$stratum), function(st) {
  s  <- per_rep[per_rep$stratum == st, ]
  ci <- boot_one(s)
  data.frame(
    stratum  = st,
    n_rep    = nrow(s),
    metric   = metrics,
    observed = vapply(metrics, function(m) mean(s[[m]], na.rm = TRUE), numeric(1)),
    ci_lo    = ci$ci_lo,
    ci_hi    = ci$ci_hi,
    ## the plain spread across replicates, which the submitted figure reports
    ## alongside the bootstrap interval -- they answer different questions.
    ## all-NA is a real case, not an error: every tandem duplication row is an
    ## upper bound, so that stratum has no point estimates to spread. min() over
    ## nothing returns Inf, which would print as a value.
    rep_min  = vapply(metrics, function(m)
                 if (all(is.na(s[[m]]))) NA_real_ else min(s[[m]], na.rm = TRUE), numeric(1)),
    rep_max  = vapply(metrics, function(m)
                 if (all(is.na(s[[m]]))) NA_real_ else max(s[[m]], na.rm = TRUE), numeric(1)),
    row.names = NULL
  )
}))

## ---- report -----------------------------------------------------------------

cat(sprintf("\nB = %d bootstrap draws, seed %d, resampling %d replicates with replacement\n\n",
            B, SEED, n_rep))

fmt <- function(m, v) {
  if (is.na(v)) return("     --")
  if (m == "mean_abs_err") sprintf("%.4f", v) else sprintf("%.1f%%", v)
}

## Every stratum present in the data is printed, in a stated order with anything
## unlisted appended. A hardcoded list silently dropped "tandem duplication"
## while it was spelled "tandem_duplication" here, and a stratum that vanishes
## from a report reads as a stratum that does not exist.
ord  <- c("ALL", "e1", "e2", "e4", "deletion", "inversion", "tandem duplication")
strata <- c(intersect(ord, res$stratum), setdiff(unique(res$stratum), ord))

for (st in strata) {
  r <- res[res$stratum == st, ]
  if (!nrow(r)) next
  cat(sprintf("%-20s\n", st))
  for (i in seq_len(nrow(r))) {
    cat(sprintf("  %-13s %8s   95%% CI [%s, %s]   replicates %s-%s\n",
                r$metric[i], fmt(r$metric[i], r$observed[i]),
                fmt(r$metric[i], r$ci_lo[i]), fmt(r$metric[i], r$ci_hi[i]),
                fmt(r$metric[i], r$rep_min[i]), fmt(r$metric[i], r$rep_max[i])))
  }
  if (all(is.na(r$observed[r$metric != "detection"])))
    cat("      (no point estimates: every row in this stratum is an upper bound)\n")
}

nb <- per_rep[per_rep$stratum == "ALL", "n_bounded", drop = TRUE]
cat(sprintf("\nbounded (tandem duplication, r not identifiable): %.0f per replicate, %s\n",
            mean(nb), "included in the error statistics above as r=2 upper bounds, not point estimates"))

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.table(res, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nwrote %s  (%d rows)\n", OUT, nrow(res)))

per_rep_out <- sub("\\.tsv$", "_per_replicate.tsv", OUT)
write.table(per_rep, per_rep_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s  (%d rows)\n", per_rep_out, nrow(per_rep)))
