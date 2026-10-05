#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})

options <- list(
  make_option(c("-f", "--file"), type = "character", help = "SV VCF file"),
  make_option(c("-p", "--purity"), type = "double", help = "Tumor purity fraction"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory"),
  make_option(c("-c", "--cnv"), type = "character", help = "FACETS directory"),
  make_option(c("-s", "--smp"), type = "character", help = "Sample name")
)
opts <- parse_args(OptionParser(option_list = options))

required <- c("file", "purity", "output", "cnv", "smp")
missing <- required[vapply(required, function(x) {
  is.null(opts[[x]]) || (is.character(opts[[x]]) && !nzchar(opts[[x]]))
}, logical(1))]
if (length(missing)) stop("Missing required options: ", paste(missing, collapse = ", "))
if (!file.exists(opts$file)) stop("SV VCF does not exist: ", opts$file)
dir.create(opts$output, recursive = TRUE, showWarnings = FALSE)

sv <- read.delim(
  opts$file,
  header = FALSE,
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE,
  check.names = FALSE
)
if (ncol(sv) != 11L) {
  stop("Expected exactly 11 VCF columns (two samples); found ", ncol(sv))
}
names(sv) <- c(
  "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
  "FORMAT", "normal", "tumor"
)

extract_info <- function(values, key) {
  matches <- regexec(paste0("(?:^|;)", key, "=([^;]+)"), values, perl = TRUE)
  parts <- regmatches(values, matches)
  vapply(parts, function(x) if (length(x) >= 2L) x[2] else NA_character_, character(1))
}

sv_class <- extract_info(sv$INFO, "SVTYPE")
is_bnd <- sv_class == "BND" | grepl("[\\[\\]]", sv$ALT, perl = TRUE)
mate_match <- regexec("[\\[\\]]([^:\\[\\]]+):(\\d+)[\\[\\]]", sv$ALT, perl = TRUE)
mate_parts <- regmatches(sv$ALT, mate_match)
mate_chr <- vapply(mate_parts, function(x) {
  if (length(x) >= 3L) x[2] else NA_character_
}, character(1))
mate_pos <- vapply(mate_parts, function(x) {
  if (length(x) >= 3L) x[3] else NA_character_
}, character(1))

pos1 <- suppressWarnings(as.integer(sv$POS))
chr2 <- ifelse(is_bnd, mate_chr, sv$CHROM)
pos2 <- suppressWarnings(as.integer(ifelse(is_bnd, mate_pos, extract_info(sv$INFO, "END"))))
keep <- sv$FILTER == "PASS" & sv_class != "INS"
if (anyNA(sv_class[keep]) || anyNA(pos1[keep]) || anyNA(chr2[keep]) || anyNA(pos2[keep])) {
  bad <- sv$ID[keep][
    is.na(sv_class[keep]) | is.na(pos1[keep]) | is.na(chr2[keep]) | is.na(pos2[keep])
  ]
  stop("Could not parse retained SV coordinates for: ", paste(bad, collapse = ", "))
}

simple <- data.frame(
  chr1 = sv$CHROM[keep],
  pos1 = pos1[keep],
  chr2 = chr2[keep],
  pos2 = pos2[keep],
  stringsAsFactors = FALSE
)
invert <- sv_class[keep] == "INV"
if (any(invert)) {
  original_pos1 <- simple$pos1[invert]
  simple$pos1[invert] <- simple$pos2[invert]
  simple$pos2[invert] <- original_pos1
}

sample_name <- opts$smp
write.table(
  simple,
  file.path(opts$output, paste0("svc_", sample_name, "_simple.txt")),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  quote = FALSE
)

fit_path <- file.path(opts$cnv, paste0(sample_name, ".RData"))
facet_path <- file.path(opts$cnv, paste0(sample_name, ".bed"))
if (!file.exists(fit_path)) stop("Missing FACETS fit: ", fit_path)
if (!file.exists(facet_path)) stop("Missing FACETS segments: ", facet_path)
fit_env <- new.env(parent = emptyenv())
load(fit_path, envir = fit_env)
if (!exists("fit", envir = fit_env, inherits = FALSE)) stop("FACETS RData lacks object 'fit': ", fit_path)
fit <- get("fit", envir = fit_env, inherits = FALSE)
if (is.null(fit$ploidy) || length(fit$ploidy) != 1L) stop("FACETS fit lacks scalar ploidy")

purity <- as.numeric(opts$purity)
if (!is.finite(purity) || purity < 0 || purity > 1) stop("Purity must be a fraction in [0, 1]")
spp <- data.frame(sample = sample_name, purity = purity, ploidy = fit$ploidy)
write.table(
  spp,
  file.path(opts$output, paste0("spp_", sample_name, ".txt")),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  quote = FALSE
)

facet <- read.delim(facet_path)
required_facet <- c("chrom", "start", "end", "tcn.em", "lcn.em")
missing_facet <- setdiff(required_facet, names(facet))
if (length(missing_facet)) {
  stop("FACETS segments lack columns: ", paste(missing_facet, collapse = ", "))
}
cnv <- facet %>%
  mutate(
    chrom = paste0("chr", chrom),
    n_major = 2,
    n_minor = 1,
    t_total = tcn.em,
    t_minor = ifelse(is.na(lcn.em), 0L, lcn.em)
  ) %>%
  select(chrom, start, end, n_major, n_minor, t_total, t_minor)
write.table(
  cnv,
  file.path(opts$output, "cnv.txt"),
  row.names = TRUE,
  col.names = FALSE,
  sep = ",",
  quote = FALSE
)

message("Wrote ", nrow(simple), " retained SVs for ", sample_name)
