#!/usr/bin/env Rscript

# Build the condition-level ledger for the VISOR three-arm comparison.
#
# This script deliberately records unavailable fair-arm estimates instead of
# dropping them. The output is the gate for per-SV scoring: every one of the
# 2,250 designed conditions must have SVCFit inputs and an assisted SVclone
# result, while the fair arm must have either a readable result file or an
# explicit failed_no_purity status (never both).

args <- commandArgs(trailingOnly = TRUE)

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: 06_build_three_arm_manifest.R \\\n  --source-root DIR --assisted-root DIR --fair-root DIR --output-dir DIR\n",
    "\n",
    "Writes the condition manifest plus overall and purity-stratified summaries.\n"
  ), file = if (status == 0L) stdout() else stderr())
  quit(status = status)
}

parse_args <- function(x) {
  if (length(x) == 1L && x %in% c("-h", "--help")) usage(0L)
  if (length(x) %% 2L != 0L) usage(2L)
  out <- list()
  for (i in seq(1L, length(x), by = 2L)) {
    key <- sub("^--", "", x[[i]])
    if (identical(key, x[[i]]) || !nzchar(key)) usage(2L)
    out[[key]] <- x[[i + 1L]]
  }
  out
}

opt <- parse_args(args)
required <- c("source-root", "assisted-root", "fair-root", "output-dir")
missing_opt <- required[!vapply(required, function(x) {
  !is.null(opt[[x]]) && nzchar(opt[[x]])
}, logical(1))]
if (length(missing_opt)) {
  stop("Missing required options: ", paste(missing_opt, collapse = ", "), call. = FALSE)
}

roots <- lapply(opt[required[1:3]], normalizePath, mustWork = TRUE)
names(roots) <- c("source", "assisted", "fair")
output_dir <- normalizePath(opt[["output-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir)

reps <- paste0("rep", seq_len(30L))
experiments <- paste0("exp", seq_len(5L))
truth_purities <- c(10L, 20L, 40L, 60L, 80L)
mixtures <- c(10L, 30L, 50L)

design <- expand.grid(
  mixture_pct = mixtures,
  truth_purity_pct = truth_purities,
  experiment = experiments,
  replicate = reps,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
design$condition <- sprintf("c50p%dm%d", design$truth_purity_pct, design$mixture_pct)
design$task_id <- seq.int(0L, nrow(design) - 1L)
design <- design[c(
  "task_id", "replicate", "experiment", "condition",
  "truth_purity_pct", "mixture_pct"
)]

nonempty <- function(path) {
  info <- file.info(path)
  file.exists(path) & !is.na(info$size) & info$size > 0
}

read_fair_status <- function(path) {
  empty <- list(status = NA_character_, purity = NA_real_, note = NA_character_)
  if (!nonempty(path)) return(empty)
  x <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) != 1L || !"facets_status" %in% names(x)) return(empty)
  list(
    status = as.character(x$facets_status[[1L]]),
    purity = if ("purity_used" %in% names(x)) suppressWarnings(as.numeric(x$purity_used[[1L]])) else NA_real_,
    note = if ("note" %in% names(x)) as.character(x$note[[1L]]) else NA_character_
  )
}

n <- nrow(design)
manifest <- design
manifest$source_dir <- manifest$assisted_dir <- manifest$fair_dir <- character(n)
manifest$svcfit_inputs_ready <- logical(n)
manifest$assisted_rdata <- manifest$fair_status_file <- manifest$fair_rdata <- character(n)
manifest$assisted_ready <- manifest$fair_rdata_ready <- logical(n)
manifest$fair_status <- manifest$fair_note <- character(n)
manifest$fair_purity_used <- rep(NA_real_, n)

for (i in seq_len(n)) {
  rep <- manifest$replicate[[i]]
  exp <- manifest$experiment[[i]]
  cond <- manifest$condition[[i]]
  rel <- file.path(rep, exp, cond)

  source_dir <- file.path(roots$source, rel)
  assisted_dir <- file.path(roots$assisted, rel)
  fair_dir <- file.path(roots$fair, rel)

  svcfit_inputs <- c(
    file.path(source_dir, "SNP", paste0("het_near_sv_", cond, ".vcf")),
    file.path(source_dir, "SNP", paste0("het_on_sv_", cond, ".vcf")),
    file.path(assisted_dir, "svtyp", paste0("svt_", cond, ".vcf")),
    file.path(source_dir, "facet", paste0(cond, ".bed"))
  )
  assisted_rdata <- file.path(
    assisted_dir, "svclone", cond, "ccube_out",
    paste0(cond, "_ccube_sv_results.RData")
  )
  fair_status_file <- file.path(
    fair_dir, "svclone_fair", paste0("facets_status_", cond, ".txt")
  )
  fair_rdata <- file.path(
    fair_dir, "svclone_fair", cond, "ccube_out",
    paste0(cond, "_ccube_sv_results.RData")
  )
  status <- read_fair_status(fair_status_file)

  manifest$source_dir[[i]] <- source_dir
  manifest$assisted_dir[[i]] <- assisted_dir
  manifest$fair_dir[[i]] <- fair_dir
  manifest$svcfit_inputs_ready[[i]] <- all(nonempty(svcfit_inputs))
  manifest$assisted_rdata[[i]] <- assisted_rdata
  manifest$assisted_ready[[i]] <- nonempty(assisted_rdata)
  manifest$fair_status_file[[i]] <- fair_status_file
  manifest$fair_status[[i]] <- status$status
  manifest$fair_purity_used[[i]] <- status$purity
  manifest$fair_note[[i]] <- status$note
  manifest$fair_rdata[[i]] <- fair_rdata
  manifest$fair_rdata_ready[[i]] <- nonempty(fair_rdata)
}

manifest$fair_outcome <- ifelse(
  manifest$fair_status == "ok" & manifest$fair_rdata_ready,
  "ready",
  ifelse(
    manifest$fair_status == "failed_no_purity" & !manifest$fair_rdata_ready,
    "failed_no_purity",
    "invalid"
  )
)
manifest$svcfit_assisted_ready <- manifest$svcfit_inputs_ready & manifest$assisted_ready
manifest$svcfit_fair_ready <- manifest$svcfit_inputs_ready & manifest$fair_outcome == "ready"
manifest$three_way_ready <- manifest$svcfit_assisted_ready & manifest$fair_outcome == "ready"

summary <- data.frame(
  metric = c(
    "designed_conditions", "unique_condition_keys", "svcfit_inputs_ready",
    "assisted_ready", "fair_ready", "fair_failed_no_purity", "fair_invalid",
    "svcfit_assisted_ready", "svcfit_fair_ready", "three_way_ready"
  ),
  value = c(
    nrow(manifest),
    length(unique(paste(manifest$replicate, manifest$experiment, manifest$condition, sep = "/"))),
    sum(manifest$svcfit_inputs_ready),
    sum(manifest$assisted_ready),
    sum(manifest$fair_outcome == "ready"),
    sum(manifest$fair_outcome == "failed_no_purity"),
    sum(manifest$fair_outcome == "invalid"),
    sum(manifest$svcfit_assisted_ready),
    sum(manifest$svcfit_fair_ready),
    sum(manifest$three_way_ready)
  ),
  stringsAsFactors = FALSE
)

purity_summary <- do.call(rbind, lapply(truth_purities, function(purity) {
  x <- manifest[manifest$truth_purity_pct == purity, , drop = FALSE]
  data.frame(
    truth_purity_pct = purity,
    designed_conditions = nrow(x),
    fair_ready = sum(x$fair_outcome == "ready"),
    fair_failed_no_purity = sum(x$fair_outcome == "failed_no_purity"),
    fair_invalid = sum(x$fair_outcome == "invalid"),
    stringsAsFactors = FALSE
  )
}))

failures <- character()
expect_equal <- function(label, actual, expected) {
  if (!identical(as.integer(actual), as.integer(expected))) {
    failures <<- c(failures, sprintf("%s: expected %d, observed %d", label, expected, actual))
  }
}
expect_equal("designed conditions", nrow(manifest), 2250L)
expect_equal("unique condition keys", summary$value[summary$metric == "unique_condition_keys"], 2250L)
expect_equal("SVCFit-ready conditions", sum(manifest$svcfit_inputs_ready), 2250L)
expect_equal("assisted-ready conditions", sum(manifest$assisted_ready), 2250L)
expect_equal("fair-ready conditions", sum(manifest$fair_outcome == "ready"), 1727L)
expect_equal("fair no-purity conditions", sum(manifest$fair_outcome == "failed_no_purity"), 523L)
expect_equal("invalid fair outcomes", sum(manifest$fair_outcome == "invalid"), 0L)
expected_fair_by_purity <- c(`10` = 12L, `20` = 365L, `40` = 450L, `60` = 450L, `80` = 450L)
for (purity in names(expected_fair_by_purity)) {
  expect_equal(
    paste0("fair-ready conditions at truth purity ", purity, "%"),
    purity_summary$fair_ready[purity_summary$truth_purity_pct == as.integer(purity)],
    expected_fair_by_purity[[purity]]
  )
}

manifest_path <- file.path(output_dir, "three_arm_condition_manifest.tsv")
summary_path <- file.path(output_dir, "three_arm_manifest_summary.tsv")
purity_summary_path <- file.path(output_dir, "three_arm_fair_availability_by_truth_purity.tsv")
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(purity_summary, purity_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

print(summary, row.names = FALSE)
cat("\nFair-arm availability by truth purity:\n")
print(purity_summary, row.names = FALSE)
cat("Manifest:", manifest_path, "\n")
cat("Summary: ", summary_path, "\n")
cat("Purity:  ", purity_summary_path, "\n")

if (length(failures)) {
  stop("Three-arm manifest gate failed:\n- ", paste(failures, collapse = "\n- "), call. = FALSE)
}
cat("Three-arm manifest gate: PASS\n")
