#!/usr/bin/env Rscript

# Repository-owned FACETS command-line driver.
#
# Contract retained from the historical project driver:
#   Rscript facet.R -s SAMPLE -p WORK_DIRECTORY
#
# It reads WORK_DIRECTORY/SAMPLE.csv.gz and writes SAMPLE.bed plus SAMPLE.RData.
# Keeping this small adapter in the workflow repository removes an undocumented
# dependency on a per-user ~/facet/facet.R file.

args <- commandArgs(trailingOnly = TRUE)

usage <- function(status = 0L) {
  cat("Usage: facet.R -s SAMPLE -p WORK_DIRECTORY\n",
      "Reads:  WORK_DIRECTORY/SAMPLE.csv.gz\n",
      "Writes: WORK_DIRECTORY/SAMPLE.bed and SAMPLE.RData\n", sep = "")
  quit(status = status)
}

if (any(args %in% c("-h", "--help"))) usage(0L)

value_after <- function(short, long) {
  pos <- which(args %in% c(short, long))
  if (length(pos) != 1L || pos == length(args)) return(NULL)
  args[[pos + 1L]]
}

sample_name <- value_after("-s", "--sample")
work_dir <- value_after("-p", "--path")
if (is.null(sample_name) || !nzchar(sample_name) ||
    is.null(work_dir) || !nzchar(work_dir)) usage(2L)

if (!requireNamespace("facets", quietly = TRUE)) {
  stop("R package 'facets' is missing; run tools/setup_facets.sh")
}

input <- file.path(work_dir, paste0(sample_name, ".csv.gz"))
bed <- file.path(work_dir, paste0(sample_name, ".bed"))
rdata <- file.path(work_dir, paste0(sample_name, ".RData"))
if (!file.exists(input)) stop("FACETS count matrix not found: ", input)

set.seed(1234)
rcmat <- facets::readSnpMatrix(input)
xx <- facets::preProcSample(rcmat, ndepth = 20, gbuild = "hg38", cval = 50)

# Preserve the project's validated retry policy for degenerate profiles.
oo <- tryCatch(
  facets::procSample(xx, cval = 150, dipLogR = NULL),
  error = function(first_error) {
    message("procSample(cval=150) failed: ", conditionMessage(first_error))
    tryCatch(
      facets::procSample(xx, cval = 300, dipLogR = NULL),
      error = function(second_error) {
        message("procSample(cval=300) failed: ", conditionMessage(second_error))
        message("Retrying with cval=150 and dipLogR=0")
        facets::procSample(xx, cval = 150, dipLogR = 0)
      }
    )
  }
)

ooo <- facets::procSample(xx, cval = 500, dipLogR = oo$dipLogR)
fit <- facets::emcncf(ooo)
dat <- fit$cncf
if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0L) {
  stop("FACETS returned no cncf segments for ", sample_name)
}

write.table(dat, bed, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE)
save(rcmat, xx, oo, ooo, fit, file = rdata)

stopifnot(file.info(bed)$size > 0, file.info(rdata)$size > 0)
message("Wrote ", bed)
message("Wrote ", rdata)
