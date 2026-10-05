#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0L) {
  cat("Usage: 08_extract_svcfit_condition.R --manifest FILE --task-id N --output-root DIR --svcfit-package DIR --truth-dir DIR --rlib DIR\n",
      file = if (status == 0L) stdout() else stderr())
  quit(status = status)
}
parse_args <- function(x) {
  if (length(x) == 1L && x %in% c("-h", "--help")) usage()
  if (length(x) %% 2L) usage(2L)
  out <- list()
  for (i in seq(1L, length(x), 2L)) out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
  out
}
opt <- parse_args(args)
required <- c("manifest", "task-id", "output-root", "svcfit-package", "truth-dir", "rlib")
missing <- required[!vapply(required, function(x) !is.null(opt[[x]]) && nzchar(opt[[x]]), logical(1))]
if (length(missing)) stop("Missing options: ", paste(missing, collapse = ", "), call. = FALSE)

paths <- lapply(opt[c("manifest", "svcfit-package", "truth-dir", "rlib")], normalizePath, mustWork = TRUE)
task_id <- suppressWarnings(as.integer(opt[["task-id"]]))
if (!is.finite(task_id) || task_id < 0L || task_id > 2249L) stop("Invalid task ID", call. = FALSE)
manifest <- read.delim(paths$manifest, stringsAsFactors = FALSE, check.names = FALSE)
row <- manifest[manifest$task_id == task_id, , drop = FALSE]
if (nrow(row) != 1L) stop("Expected one manifest row for task ", task_id, call. = FALSE)

out_dir <- file.path(opt[["output-root"]], "svcfit_events", row$replicate, row$experiment)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(out_dir, paste0(row$condition, ".tsv"))
components_tsv <- file.path(out_dir, paste0(row$condition, ".components.tsv"))
status_tsv <- file.path(out_dir, paste0(row$condition, ".status.tsv"))
if (file.exists(out_tsv) && file.info(out_tsv)$size > 0L && Sys.getenv("FORCE_SVCFIT_EXTRACT", "0") != "1") {
  cat("Already complete:", out_tsv, "\n")
  quit(status = 0L)
}

.libPaths(c(paths$rlib, .libPaths()))
if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools unavailable in --rlib", call. = FALSE)
devtools::load_all(paths[["svcfit-package"]], quiet = TRUE)

cond <- row$condition
source_dir <- row$source_dir
assisted_dir <- row$assisted_dir
inputs <- list(
  p_het = file.path(source_dir, "SNP", paste0("het_near_sv_", cond, ".vcf")),
  p_onsv = file.path(source_dir, "SNP", paste0("het_on_sv_", cond, ".vcf")),
  p_sv = file.path(assisted_dir, "svtyp", paste0("svt_", cond, ".vcf")),
  p_cnv = file.path(source_dir, "facet", paste0(cond, ".bed"))
)
if (any(!file.exists(unlist(inputs)))) stop("Missing SVCFit input for task ", task_id, call. = FALSE)

fit <- do.call(run_svcfit, c(inputs, list(
  samp = cond, exper = row$experiment, QUAL_thresh = 100,
  min_alt = 2, tum_only = FALSE
)))
truth <- load_truth(paths[["truth-dir"]], overlap = row$experiment != "exp1")
final <- attach_truth(fit$svcf, truth)
if (!nrow(final)) stop("SVCFit produced no truth-annotated events", call. = FALSE)

is_bnd <- final$classification == "BND" | grepl("BND", final$ID)
second_pos <- ifelse(is_bnd, final$pos2, final$END)
mutation_id <- paste(final$CHROM, final$POS, final$chr2, second_pos, sep = ":")
if (anyNA(mutation_id) || any(!nzchar(mutation_id))) stop("Malformed SVCFit event key", call. = FALSE)
purity <- as.numeric(row$truth_purity_pct) / 100
mixture <- as.numeric(row$mixture_pct) / 100
true_svcf <- ifelse(final$clone == "clonal", purity,
                    ifelse(final$clone == "2", purity * mixture,
                           ifelse(final$clone == "sub", purity * (1 - mixture), NA_real_)))
components <- data.frame(
  task_id = task_id, replicate = row$replicate, experiment = row$experiment,
  condition = cond, truth_purity_pct = row$truth_purity_pct,
  mixture_pct = row$mixture_pct, mutation_id = mutation_id,
  svcfit_svcf = as.numeric(final$final_svcf), true_svcf = true_svcf,
  signed_error = true_svcf - as.numeric(final$final_svcf),
  abs_error = abs(true_svcf - as.numeric(final$final_svcf)),
  classification = as.character(final$classification),
  zygosity = as.character(final$zygosity), clone = as.character(final$clone),
  svcfit_status = as.character(final$svcf_status),
  stringsAsFactors = FALSE
)
components$score_ready <- is.finite(components$svcfit_svcf) & is.finite(components$true_svcf)

# SVCFit can retain paired components for a single symbolic event (notably
# INV). Preserve those rows, then average them once per strand-free event key,
# matching the event-level rule used by 07_extract_svclone_arms.R.
groups <- split(seq_len(nrow(components)), components$mutation_id)
events <- do.call(rbind, lapply(groups, function(idx) {
  x <- components[idx, , drop = FALSE]
  invariant <- c("classification", "zygosity", "clone", "true_svcf")
  bad <- invariant[vapply(invariant, function(nm) length(unique(x[[nm]])) > 1L, logical(1))]
  if (length(bad)) stop("Inconsistent paired-component metadata for ", x$mutation_id[[1L]],
                        ": ", paste(bad, collapse = ", "), call. = FALSE)
  estimate <- if (all(is.na(x$svcfit_svcf))) NA_real_ else mean(x$svcfit_svcf, na.rm = TRUE)
  data.frame(
    task_id = x$task_id[[1L]], replicate = x$replicate[[1L]],
    experiment = x$experiment[[1L]], condition = x$condition[[1L]],
    truth_purity_pct = x$truth_purity_pct[[1L]], mixture_pct = x$mixture_pct[[1L]],
    mutation_id = x$mutation_id[[1L]], component_count = nrow(x),
    svcfit_svcf = estimate, true_svcf = x$true_svcf[[1L]],
    signed_error = x$true_svcf[[1L]] - estimate,
    abs_error = abs(x$true_svcf[[1L]] - estimate),
    classification = x$classification[[1L]], zygosity = x$zygosity[[1L]],
    clone = x$clone[[1L]],
    svcfit_status = paste(unique(x$svcfit_status), collapse = ";"),
    score_ready = is.finite(estimate) && is.finite(x$true_svcf[[1L]]),
    stringsAsFactors = FALSE
  )
}))
rownames(events) <- NULL
if (anyDuplicated(events$mutation_id)) stop("Event aggregation left duplicate keys", call. = FALSE)

tmp <- tempfile(pattern = paste0(cond, "."), tmpdir = out_dir)
write.table(events, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (!file.rename(tmp, out_tsv)) stop("Atomic output rename failed", call. = FALSE)
tmp_components <- tempfile(pattern = paste0(cond, ".components."), tmpdir = out_dir)
write.table(components, tmp_components, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (!file.rename(tmp_components, components_tsv)) stop("Atomic component-output rename failed", call. = FALSE)
status <- data.frame(task_id = task_id, components = nrow(components), events = nrow(events),
                     score_ready = sum(events$score_ready), status = "ok")
write.table(status, status_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Done: SVCFit extraction [%s/%s/%s] events=%d score_ready=%d\n",
            row$replicate, row$experiment, cond, nrow(events), sum(events$score_ready)))
