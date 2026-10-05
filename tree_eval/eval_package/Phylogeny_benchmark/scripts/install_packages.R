#!/usr/bin/env Rscript
#
# Install the R dependencies for this evaluation package.
# Run once on the local machine:   Rscript scripts/install_packages.R
#
# Extends the 3-package helper from tree_eval_package/a2_coverage_correlation/
# to the full set evaluate_downstream.R actually loads.

pkgs <- c("optparse", "dplyr", "readr", "tidyr", "ggplot2", "patchwork", "lme4")

repo <- "https://cloud.r-project.org"
to_install <- setdiff(pkgs, rownames(installed.packages()))

if (length(to_install) > 0L) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = repo)
} else {
  message("All required packages already installed.")
}

ok_all <- TRUE
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  ok_all <- ok_all && ok
  message(sprintf("  %-10s %s", p, if (ok) "OK" else "MISSING"))
}

if (!ok_all) quit(status = 1L)
