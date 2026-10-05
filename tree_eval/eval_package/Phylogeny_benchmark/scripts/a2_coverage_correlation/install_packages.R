#!/usr/bin/env Rscript
#
# Install R package dependencies for the A2 coverage correlation analysis.
# Run once on the cluster:   Rscript install_packages.R

pkgs <- c("dplyr", "ggplot2", "lme4")

repo <- "https://cloud.r-project.org"
to_install <- setdiff(pkgs, rownames(installed.packages()))

if (length(to_install) > 0L) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = repo)
} else {
  message("All required packages already installed.")
}

for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  message(sprintf("  %-10s %s", p, if (ok) "OK" else "MISSING"))
}
