#!/usr/bin/env Rscript
## Which form should a HEMIZYGOUS DELETION use?
##
## Scores the three candidates against known cellular fraction on every detected
## deletion in all 45 conditions, so the choice rests on the whole purity x
## mixture grid rather than one point:
##
##   h1     = cn_bar*VAF - (cn_bar - 1)   SV before CNV   (what the sign rule picks)
##   h2     = cn_bar*VAF                  CNV before SV
##   depth  = 1 - cn_bar                  carrier cells contribute no alleles
##
## Reads the scored table written by 13_chrx_score_svcfit.R. Routing-independent:
## it recomputes all three from cn_bar and VAF, so the dispatch under test cannot
## bias the comparison.

suppressMessages({ library(dplyr) })

CHRX_DIR <- Sys.getenv("CHRX_DIR", unset = ".")
## Coverage-named, as 13 writes it. COV is required; see the note in chrx_common.sh.
.cov <- Sys.getenv("COV", unset = "")
if (!nzchar(.cov)) stop("COV is not set, and there is no default. Run e.g. COV=50 ...",
                        call. = FALSE)
ctag <- paste0("c", .cov)
## Replicate tag, as 13 writes it; empty at REP=0. See chrx_common.sh.
REP  <- Sys.getenv("REP", unset = "0")
rtag <- if (REP == "0") "" else paste0("_rep", REP)
f <- file.path(CHRX_DIR, paste0("scoring", rtag),
               paste0("chrx_svcfit_scores_", ctag, ".tsv"))
d <- read.delim(f, stringsAsFactors = FALSE)

dl <- d %>%
  filter(sv_class == "deletion", detected, is.finite(cn_bar), is.finite(vaf)) %>%
  mutate(h1    = cn_bar * vaf - (cn_bar - 1),
         h2    = cn_bar * vaf,
         depth = 1 - cn_bar)

## Long form, so every summary below scores the three identically.
long <- bind_rows(
  dl %>% mutate(form = "h1",    est = h1),
  dl %>% mutate(form = "h2",    est = h2),
  dl %>% mutate(form = "depth", est = depth)
) %>% mutate(err = est - cellular_fraction)

summ <- function(g) g %>% summarise(
  n = n(),
  med_abs = median(abs(err)),
  mean_abs = mean(abs(err)),
  within_05 = 100 * mean(abs(err) <= 0.05),
  ## an estimate outside (0,1] is not a cellular fraction at all
  infeasible = 100 * mean(est <= 0 | est > 1),
  .groups = "drop")

cat("=== deletions, all 45 conditions ===\n")
cat("total detected deletions:", nrow(dl), "\n\n")

cat("by experiment x form:\n")
print(as.data.frame(summ(long %>% group_by(experiment, form))), row.names = FALSE)

cat("\nby purity x form (all experiments):\n")
print(as.data.frame(summ(long %>% group_by(purity_pct, form))), row.names = FALSE)

cat("\nby mixture x form (all experiments):\n")
print(as.data.frame(summ(long %>% group_by(mixture_pct, form))), row.names = FALSE)

cat("\nby experiment x purity, best form by median |err|:\n")
best <- long %>% group_by(experiment, purity_pct, form) %>%
  summarise(med = median(abs(err)), .groups = "drop") %>%
  group_by(experiment, purity_pct) %>%
  slice_min(med, n = 1, with_ties = FALSE) %>% ungroup()
print(as.data.frame(best), row.names = FALSE)

## Is the winner stable across the grid, or an artefact of one corner?
cat("\nwinner count across the 45 experiment x condition cells:\n")
cells <- long %>% group_by(experiment, condition, form) %>%
  summarise(med = median(abs(err)), .groups = "drop") %>%
  group_by(experiment, condition) %>%
  slice_min(med, n = 1, with_ties = FALSE) %>% ungroup()
print(as.data.frame(table(cells$experiment, cells$form)), row.names = FALSE)

## The sign rule's own behaviour: how often is h1 even close to the zero threshold?
cat("\nsign-rule margin (h1 vs its hard-zero threshold):\n")
print(as.data.frame(dl %>% group_by(experiment) %>% summarise(
  n = n(), med_h1 = median(h1), frac_h1_pos = 100 * mean(h1 > 0),
  frac_near_zero = 100 * mean(abs(h1) < 0.05), .groups = "drop")), row.names = FALSE)
