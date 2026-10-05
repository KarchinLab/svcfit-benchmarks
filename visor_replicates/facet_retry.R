library(facets)
library(optparse)
library(readr)

option_list <- list(
  make_option(c("-s", "--sample"),
              type = "character",
              help = "sample name"),
  make_option(c("-p", "--path"),
              type = "character",
              help = "working directory")
)

parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opts   <- parse_args(parser)

set.seed(1234)
SNP_file <- paste0(opts$path, "/", opts$sample, ".csv.gz")
rcmat <- readSnpMatrix(SNP_file)
xx <- preProcSample(rcmat, ndepth = 20, gbuild = "hg38", cval = 50)

# First attempt: standard call with auto dipLogR detection via clustersegs.
# This can crash with "non-numeric argument to mathematical function" when the
# copy-number profile produces a degenerate cluster (lorclust[[jhi]]$valor is NULL).
oo <- tryCatch(
  procSample(xx, cval = 150, dipLogR = NULL),
  error = function(e) {
    message("[facet_retry] procSample(cval=150, dipLogR=NULL) failed: ", conditionMessage(e))
    message("[facet_retry] Retrying with cval=300 ...")
    tryCatch(
      procSample(xx, cval = 300, dipLogR = NULL),
      error = function(e2) {
        message("[facet_retry] procSample(cval=300, dipLogR=NULL) failed: ", conditionMessage(e2))
        message("[facet_retry] Falling back to dipLogR=0 (safe default for diploid regions).")
        procSample(xx, cval = 150, dipLogR = 0)
      }
    )
  }
)

ooo <- procSample(xx, cval = 500, dipLogR = oo$dipLogR)
fit <- emcncf(ooo)
dat <- fit$cncf
write_delim(dat, paste0(opts$path, "/", opts$sample, ".bed"),
            quote = "none", col_names = TRUE, delim = "\t")
save(rcmat, xx, oo, ooo, fit, file = paste0(opts$path, "/", opts$sample, ".RData"))
