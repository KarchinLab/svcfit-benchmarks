#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# evaluate_downstream.R — CLI wrapper around evaluate_downstream.Rmd
#
# Cluster size pools pre + post SV counts
# into a single per-replicate total before computing that correlation, matching
# the Supplement's stated methods ("cluster size is computed as the pair-total
# count, combining pre- and post-treatment timepoints"). It reproduces the
# the Supplement's values exactly: rho = -0.5074 / -0.4021. See
# ../docs/WHY_evaluate_downstream_Rmd_is_canonical.md for the full diagnosis.
#
# To avoid the older logic being run by accident, only the .Rmd now ships. This
# wrapper preserves the command-line interface that run_all.sh,
# run_svcfit_and_evaluate.sh and resubmit_skipped.sh expect, and drives the
# canonical .Rmd by overriding its parameter chunk.
#
# Usage (unchanged):
#   Rscript evaluate_downstream.R --work_dir DIR --input_dir DIR \
#     --truth_dir DIR --out_dir DIR [--cov 50] [--n_boots 5] [--boot N] \
#     [--expected_cases 100] [--expected_no_tree 8] \
#     [--expected_correct_topologies 80] [--purities 20,40,60,80] \
#     [--experiments exp1,exp2,exp3,exp4,exp5]
# ---------------------------------------------------------------------------

suppressMessages(library(knitr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag)
  args[i + 1]
}

script_dir <- {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) normalizePath(dirname(sub("^--file=", "", fa[1]))) else getwd()
}

work_dir  <- get_arg("--work_dir")
input_dir <- get_arg("--input_dir", work_dir)
truth_dir <- get_arg("--truth_dir")
out_dir   <- get_arg("--out_dir")
if (is.null(work_dir) || is.null(truth_dir) || is.null(out_dir))
  stop("--work_dir, --truth_dir and --out_dir are all required")

cov     <- as.integer(get_arg("--cov", "50"))
n_boots <- as.integer(get_arg("--n_boots", "5"))
boot    <- get_arg("--boot", NA)
boot    <- if (is.na(boot)) NA_integer_ else as.integer(boot)
expected_cases <- as.integer(get_arg("--expected_cases", if (is.na(boot)) "100" else "20"))
expected_no_tree <- as.integer(get_arg("--expected_no_tree", NA))
expected_correct_topologies <- as.integer(get_arg("--expected_correct_topologies", NA))
purities <- as.integer(strsplit(get_arg("--purities", "20,40,60,80"), ",", fixed = TRUE)[[1]])
experiments <- strsplit(get_arg("--experiments", "exp1,exp2,exp3,exp4,exp5"), ",", fixed = TRUE)[[1]]
if (is.na(expected_cases) || expected_cases < 1L) stop("--expected_cases must be a positive integer")
if (!length(purities) || anyNA(purities) || any(!purities %in% c(20L, 40L, 60L, 80L)))
  stop("--purities must be a comma-separated subset of 20,40,60,80")
if (!length(experiments) || any(!experiments %in% paste0("exp", 1:5)))
  stop("--experiments must be a comma-separated subset of exp1,...,exp5")

rmd <- file.path(script_dir, "evaluate_downstream.Rmd")
if (!file.exists(rmd)) stop("Canonical notebook not found: ", rmd)

# Replace the notebook's hardcoded parameter chunk with the CLI values, so the
# analysis code below it is used verbatim and cannot drift from the .Rmd.
src <- readLines(rmd, warn = FALSE)
p_start <- grep("^```\\{r parameters\\}", src)[1]
p_end   <- grep("^```\\s*$", src)
p_end   <- p_end[p_end > p_start][1]
if (is.na(p_start) || is.na(p_end)) stop("Could not locate the parameters chunk in ", rmd)

parameter_lines <- c(
  sprintf('work_dir   <- %s', deparse(work_dir)),
  sprintf('input_dir  <- %s', deparse(input_dir)),
  sprintf('truth_dir  <- %s', deparse(truth_dir)),
  sprintf('out_dir    <- %s', deparse(out_dir)),
  sprintf('cov        <- %dL', cov),
  sprintf('n_boots    <- %dL', n_boots),
  sprintf('expected_cases <- %dL', expected_cases),
  sprintf('purities <- c(%s)', paste(sprintf("%dL", purities), collapse = ",")),
  sprintf('experiments <- c(%s)', paste(vapply(experiments, deparse, character(1)), collapse = ",")),
  if (is.na(expected_no_tree)) 'expected_no_tree <- NA_integer_' else sprintf('expected_no_tree <- %dL', expected_no_tree),
  if (is.na(expected_correct_topologies)) 'expected_correct_topologies <- NA_integer_' else sprintf('expected_correct_topologies <- %dL', expected_correct_topologies),
  if (is.na(boot)) 'boot       <- NA_integer_' else sprintf('boot       <- %dL', boot),
  'pre_stage  <- "pre_BAT"',
  'post_stage <- "on_BAT"',
  'dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)'
)
src <- c(src[seq_len(p_start)], parameter_lines, src[p_end:length(src)])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tmp_rmd <- file.path(tempdir(), "evaluate_downstream_cli.Rmd")
writeLines(src, tmp_rmd)

cat("Running canonical evaluate_downstream.Rmd via CLI wrapper\n")
cat("  work_dir : ", work_dir, "\n  input_dir: ", input_dir, "\n  truth_dir: ", truth_dir,
    "\n  out_dir  : ", out_dir, "\n  cov=", cov, " n_boots=", n_boots, "\n", sep = "")

old <- setwd(out_dir); on.exit(setwd(old))
knitr::opts_chunk$set(fig.path = "knitr_figs/", dev = "png", dpi = 150,
                      error = FALSE)
rendered <- knitr::knit(tmp_rmd, output = file.path(out_dir, "evaluate_downstream.md"))
if (!file.exists(rendered) || file.info(rendered)$size <= 0)
  stop("Evaluation did not create a nonempty report")
cat("Done. Outputs in ", out_dir, "\n", sep = "")
