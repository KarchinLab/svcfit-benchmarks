#!/usr/bin/env Rscript
##
## Run SVCFit over the prostate mixture grid with the hemizygous chrX treatment, and write the
## per-sample results. SVCFit only -- SVclone is not run.
##
## It evaluates the `setup` and `helpers` chunks of 06_svcfit_replicates_shared_sv.Rmd, so analyze()
## has one definition.
##
## Paths come from the environment (PM_PUB_ROOT, PM_REP_BASE); see the Rmd's setup chunk.
##
## 4m and 5m are excluded by default: they are built from chromosome-split sources covering
## autosomes only, so they contain no X and have no cn_bar tables. The same driver produces their
## autosome-only values with
##
##   PM_CONDITIONS=4m,5m PM_OUT_NAME=svcfit_45_all.tsv Rscript 06a_run_svcfit_chrx.R
##
## With no cn_bar table the hemizygous arguments are NULL, so those rows use the autosomal model only.

suppressMessages(library(dplyr))

RMD <- Sys.getenv("PM_RMD", file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                                      "06_svcfit_replicates_shared_sv.Rmd"))
if (!file.exists(RMD)) stop("Rmd not found: ", RMD, call. = FALSE)

## ---- evaluate only the setup and helpers chunks --------------------------------------------------
src <- readLines(RMD, warn = FALSE)
starts <- grep("^```\\{r", src)
ends   <- grep("^```$", src)
want   <- c("setup", "helpers")
code   <- character(0)
for (s in starts) {
  nm <- sub("^```\\{r[, ]*([A-Za-z0-9_.-]*).*$", "\\1", src[s])
  if (!nm %in% want) next
  e <- ends[ends > s][1]
  if (!is.na(e) && e > s + 1) code <- c(code, src[(s + 1):(e - 1)])
}
if (!length(code)) stop("no setup/helpers chunks found in ", RMD, call. = FALSE)
cat(sprintf("[%s] evaluating %d lines from the Rmd's setup+helpers chunks\n",
            format(Sys.time()), length(code)))
eval(parse(text = paste(code, collapse = "\n")), envir = globalenv())

for (o in c("analyze", "rep_base_dir", "exp_3"))
  if (!exists(o)) stop("the Rmd's chunks did not define ", o, call. = FALSE)

OUT <- Sys.getenv("PM_SVCFIT_OUT", file.path(rep_base_dir, "svcfit_chrx"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REPS <- as.integer(strsplit(Sys.getenv("PM_REPS", paste(1:30, collapse = ",")), ",")[[1]])

cat(sprintf("[%s] rep_base_dir: %s\n", format(Sys.time()), rep_base_dir))
CONDS <- strsplit(Sys.getenv("PM_CONDITIONS", paste(exp_3, collapse = ",")), ",")[[1]]
OUT_NAME <- Sys.getenv("PM_OUT_NAME", "svcfit_chrx_all.tsv")
cat(sprintf("[%s] conditions:   %s\n", format(Sys.time()), paste(CONDS, collapse = " ")))
cat(sprintf("[%s] output:       %s\n", format(Sys.time()), OUT))

## ---- run -----------------------------------------------------------------------------------------
all <- list(); n_ok <- 0L; n_fail <- 0L
for (r in REPS) {
  bd <- file.path(rep_base_dir, paste0("rep", r))
  if (!dir.exists(bd)) { cat(sprintf("  rep%-3d MISSING %s\n", r, bd)); next }
  for (ex in CONDS) {
    res <- analyze(ex, bd)          # analyze() traps its own errors and returns NULL
    if (is.null(res) || !nrow(res)) { n_fail <- n_fail + 1L; next }
    res$replicate <- r
    res$condition <- ex
    all[[length(all) + 1]] <- res
    n_ok <- n_ok + 1L
  }
  cat(sprintf("  rep%-3d done (%d ok, %d failed so far)\n", r, n_ok, n_fail))
}
if (!length(all)) stop("no samples produced results", call. = FALSE)

svcf <- bind_rows(all)
f <- file.path(OUT, OUT_NAME)
write.table(svcf, f, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\n[%s] wrote %s  (%d rows, %d samples ok, %d failed)\n",
            format(Sys.time()), f, nrow(svcf), n_ok, n_fail))

## ---- what actually happened on chrX ----------------------------------------------------------------
## The point of the exercise. If chrX rows are absent or all unresolved, the run did not do what it
## was configured to do, and that has to be visible here rather than inferred later from a figure.
if (!"CHROM" %in% names(svcf)) {
  cat("NOTE: no CHROM column in the result; cannot report the chrX breakdown\n")
} else {
  x <- svcf[svcf$CHROM %in% c("chrX", "X"), , drop = FALSE]
  cat(sprintf("\nchrX rows: %d of %d total\n", nrow(x), nrow(svcf)))
  if (nrow(x)) {
    if ("svcf_status" %in% names(x)) { cat("chrX svcf_status:\n"); print(table(x$svcf_status, useNA = "ifany")) }
    if ("sv_cnv_order" %in% names(x)) { cat("chrX sv_cnv_order:\n"); print(table(x$sv_cnv_order, useNA = "ifany")) }
    if ("svcf_is_bound" %in% names(x)) cat(sprintf("chrX upper bounds: %d\n", sum(x$svcf_is_bound %in% TRUE)))
  } else if (all(CONDS %in% c("4m", "5m"))) {
    cat("no chrX rows, which is correct for 4m/5m: they are built from autosome-only sources.\n")
  } else {
    cat("** NO chrX ROWS. chr_lst, hemizygous_chr or the cn_bar tables are not reaching the estimator. **\n")
  }
}
