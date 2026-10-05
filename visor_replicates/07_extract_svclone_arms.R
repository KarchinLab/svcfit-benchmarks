#!/usr/bin/env Rscript

# Extract the saved assisted and fair SVclone estimates into auditable tables.
# This stage does not score accuracy and does not use truth to alter the fair
# estimates. It validates the RData schema, preserves native side-level CCF and
# multiplicity fields, and constructs the fair-versus-assisted shared-SV set.

args <- commandArgs(trailingOnly = TRUE)

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: 07_extract_svclone_arms.R \\\n",
    "  --manifest FILE --output-dir DIR\n"
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
required_opt <- c("manifest", "output-dir")
missing_opt <- required_opt[!vapply(required_opt, function(x) {
  !is.null(opt[[x]]) && nzchar(opt[[x]])
}, logical(1))]
if (length(missing_opt)) {
  stop("Missing required options: ", paste(missing_opt, collapse = ", "), call. = FALSE)
}

manifest_path <- normalizePath(opt$manifest, mustWork = TRUE)
output_dir <- normalizePath(opt[["output-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir)

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_manifest <- c(
  "task_id", "replicate", "experiment", "condition", "truth_purity_pct",
  "mixture_pct", "assisted_rdata", "assisted_ready", "fair_rdata",
  "fair_purity_used", "fair_outcome"
)
missing_columns <- setdiff(required_manifest, names(manifest))
if (length(missing_columns)) {
  stop("Manifest is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (nrow(manifest) != 2250L || anyDuplicated(manifest$task_id)) {
  stop("Manifest must contain 2,250 unique task IDs", call. = FALSE)
}
if (!all(manifest$assisted_ready)) stop("Manifest contains unavailable assisted results", call. = FALSE)
if (sum(manifest$fair_outcome == "ready") != 1727L ||
    sum(manifest$fair_outcome == "failed_no_purity") != 523L ||
    any(!manifest$fair_outcome %in% c("ready", "failed_no_purity"))) {
  stop("Manifest fair outcomes do not satisfy the 1,727/523 gate", call. = FALSE)
}

required_ssm <- c(
  "mutation_id", "purity", "ccube_ccf1", "ccube_ccf2",
  "ccube_mult1", "ccube_mult2", "ccube_ccf_mean"
)

normalize_mutation_id <- function(x) {
  x <- sub(":[+-]_", ":", as.character(x))
  sub(":[+-]$", "", x)
}

load_arm <- function(path, meta, arm, expected_purity = NA_real_) {
  env <- new.env(parent = emptyenv())
  loaded <- tryCatch(load(path, envir = env), error = function(e) e)
  if (inherits(loaded, "error")) {
    stop("Could not load ", arm, " RData for task ", meta$task_id,
         ": ", conditionMessage(loaded), call. = FALSE)
  }
  if (!exists("doubleBreakPtsRes", envir = env, inherits = FALSE)) {
    stop("Missing doubleBreakPtsRes in ", arm, " RData for task ", meta$task_id, call. = FALSE)
  }
  result <- env$doubleBreakPtsRes
  if (!is.list(result) || is.null(result$ssm) || !is.data.frame(result$ssm)) {
    stop("Invalid doubleBreakPtsRes$ssm in ", arm, " RData for task ", meta$task_id, call. = FALSE)
  }
  ssm <- result$ssm
  missing_ssm <- setdiff(required_ssm, names(ssm))
  if (length(missing_ssm)) {
    stop("Missing ", arm, " SSM columns for task ", meta$task_id,
         ": ", paste(missing_ssm, collapse = ", "), call. = FALSE)
  }
  if (!nrow(ssm)) stop("Empty ", arm, " SSM table for task ", meta$task_id, call. = FALSE)

  filtered_path <- file.path(dirname(dirname(path)), paste0(meta$condition, "_filtered_svs.tsv"))
  filtered <- read.delim(filtered_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_filtered <- c("chr1", "pos1", "dir1", "chr2", "pos2", "dir2", "classification")
  missing_filtered <- setdiff(required_filtered, names(filtered))
  if (length(missing_filtered)) stop("Invalid filtered-SV table for task ", meta$task_id, call. = FALSE)
  filtered$mutation_id_raw <- paste0(filtered$chr1, ":", filtered$pos1, ":", filtered$dir1,
                                     "_", filtered$chr2, ":", filtered$pos2, ":", filtered$dir2)
  match_row <- match(as.character(ssm$mutation_id), filtered$mutation_id_raw)
  if (anyNA(match_row)) stop("RData events missing from filtered-SV table for task ", meta$task_id, call. = FALSE)

  purity_values <- unique(as.numeric(ssm$purity[is.finite(as.numeric(ssm$purity))]))
  if (length(purity_values) != 1L) {
    stop("Expected one finite ", arm, " purity for task ", meta$task_id, call. = FALSE)
  }
  purity_used <- purity_values[[1L]]
  # SVclone serializes purity at six decimal places in the RData while the
  # status ledger retains the full FACETS value.
  if (is.finite(expected_purity) && abs(purity_used - expected_purity) > 1e-6) {
    stop("Fair status/RData purity mismatch for task ", meta$task_id,
         ": status=", expected_purity, ", RData=", purity_used, call. = FALSE)
  }

  mutation_id <- normalize_mutation_id(ssm$mutation_id)
  if (anyNA(mutation_id) || any(!nzchar(mutation_id))) {
    stop("Missing normalized mutation ID in ", arm, " task ", meta$task_id, call. = FALSE)
  }
  ccf1 <- as.numeric(ssm$ccube_ccf1)
  ccf2 <- as.numeric(ssm$ccube_ccf2)
  ccf_native <- rowMeans(cbind(ccf1, ccf2))
  ccf_capped <- pmin(pmax(ccf_native, 0), 1)

  data.frame(
    task_id = meta$task_id,
    replicate = meta$replicate,
    experiment = meta$experiment,
    condition = meta$condition,
    truth_purity_pct = meta$truth_purity_pct,
    mixture_pct = meta$mixture_pct,
    arm = arm,
    component_index = seq_len(nrow(ssm)),
    mutation_id_raw = as.character(ssm$mutation_id),
    mutation_id = mutation_id,
    chr1 = filtered$chr1[match_row], pos1 = as.integer(filtered$pos1[match_row]),
    chr2 = filtered$chr2[match_row], pos2 = as.integer(filtered$pos2[match_row]),
    classification = filtered$classification[match_row],
    purity_used = purity_used,
    ccube_ccf1 = ccf1,
    ccube_ccf2 = ccf2,
    ccube_ccf_native = ccf_native,
    ccube_ccf_capped = ccf_capped,
    svcf_native = ccf_capped * purity_used,
    ccube_mult1 = as.numeric(ssm$ccube_mult1),
    ccube_mult2 = as.numeric(ssm$ccube_mult2),
    ccube_cluster_ccf = as.numeric(ssm$ccube_ccf_mean),
    stringsAsFactors = FALSE
  )
}

assisted <- vector("list", nrow(manifest))
fair_rows <- which(manifest$fair_outcome == "ready")
fair <- vector("list", length(fair_rows))
fair_index <- 0L

for (i in seq_len(nrow(manifest))) {
  meta <- manifest[i, , drop = FALSE]
  assisted[[i]] <- load_arm(meta$assisted_rdata, meta, "assisted")
  if (meta$fair_outcome == "ready") {
    fair_index <- fair_index + 1L
    fair[[fair_index]] <- load_arm(
      meta$fair_rdata, meta, "fair", as.numeric(meta$fair_purity_used)
    )
  }
  if (i %% 100L == 0L || i == nrow(manifest)) {
    cat(sprintf("Loaded %d/%d conditions\n", i, nrow(manifest)))
  }
}

assisted <- do.call(rbind, assisted)
fair <- do.call(rbind, fair)
components <- rbind(assisted, fair)
rownames(components) <- NULL

# SVclone expands some symbolic events (notably INV) into multiple signed
# components. The strand-free key is the event key used by SVCFit, so retain
# every native component above and form exactly one event-level estimate by
# averaging its component estimates. This avoids a many-to-many join and avoids
# giving multi-component inversions extra weight in event-level accuracy.
aggregate_events <- function(x) {
  group <- c("task_id", "mutation_id")
  metadata <- c(
    "replicate", "experiment", "condition", "truth_purity_pct",
    "mixture_pct", "arm", "chr1", "pos1", "chr2", "pos2", "classification"
  )
  numeric_values <- c(
    "purity_used", "ccube_ccf1", "ccube_ccf2", "ccube_ccf_native",
    "ccube_ccf_capped", "svcf_native", "ccube_mult1", "ccube_mult2",
    "ccube_cluster_ccf"
  )
  first <- x[!duplicated(x[group]), c(group, metadata), drop = FALSE]
  means <- aggregate(x[numeric_values], by = x[group], FUN = function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  })
  raw_ids <- aggregate(x["mutation_id_raw"], by = x[group], FUN = function(z) {
    paste(unique(z), collapse = ";")
  })
  counts <- aggregate(x["component_index"], by = x[group], FUN = length)
  names(counts)[names(counts) == "component_index"] <- "component_count"
  out <- Reduce(function(a, b) merge(a, b, by = group, sort = FALSE),
                list(first, raw_ids, counts, means))
  out[order(out$task_id), ]
}

assisted_events <- aggregate_events(assisted)
fair_events <- aggregate_events(fair)
events <- rbind(assisted_events, fair_events)
rownames(events) <- NULL

key <- c(
  "task_id", "replicate", "experiment", "condition",
  "truth_purity_pct", "mixture_pct", "mutation_id"
)
value_columns <- c(
  "mutation_id_raw", "component_count", "chr1", "pos1", "chr2", "pos2", "classification",
  "purity_used", "ccube_ccf1", "ccube_ccf2",
  "ccube_ccf_native", "ccube_ccf_capped", "svcf_native",
  "ccube_mult1", "ccube_mult2", "ccube_cluster_ccf"
)
pairwise <- merge(
  assisted_events[c(key, value_columns)],
  fair_events[c(key, value_columns)],
  by = key,
  suffixes = c("_assisted", "_fair"),
  all = FALSE,
  sort = FALSE
)
pairwise$delta_svcf_fair_minus_assisted <-
  pairwise$svcf_native_fair - pairwise$svcf_native_assisted

count_by_task <- function(x, name) {
  tab <- table(x$task_id)
  out <- data.frame(task_id = as.integer(names(tab)), value = as.integer(tab))
  names(out)[[2L]] <- name
  out
}
condition_summary <- manifest[c(
  "task_id", "replicate", "experiment", "condition", "truth_purity_pct",
  "mixture_pct", "fair_outcome", "fair_purity_used"
)]
condition_summary <- merge(
  condition_summary, count_by_task(assisted_events, "assisted_sv_count"),
  by = "task_id", all.x = TRUE, sort = FALSE
)
condition_summary <- merge(
  condition_summary, count_by_task(fair_events, "fair_sv_count"),
  by = "task_id", all.x = TRUE, sort = FALSE
)
condition_summary <- merge(
  condition_summary, count_by_task(pairwise, "shared_sv_count"),
  by = "task_id", all.x = TRUE, sort = FALSE
)
for (column in c("assisted_sv_count", "fair_sv_count", "shared_sv_count")) {
  condition_summary[[column]][is.na(condition_summary[[column]])] <- 0L
}
condition_summary$shared_fraction_assisted <- with(
  condition_summary,
  ifelse(assisted_sv_count > 0L, shared_sv_count / assisted_sv_count, NA_real_)
)
condition_summary$shared_fraction_fair <- with(
  condition_summary,
  ifelse(fair_sv_count > 0L, shared_sv_count / fair_sv_count, NA_real_)
)
condition_summary <- condition_summary[order(condition_summary$task_id), ]

if (nrow(condition_summary) != 2250L ||
    any(condition_summary$assisted_sv_count <= 0L) ||
    any(condition_summary$fair_outcome == "ready" & condition_summary$fair_sv_count <= 0L) ||
    any(condition_summary$fair_outcome == "failed_no_purity" & condition_summary$fair_sv_count != 0L)) {
  stop("Extracted condition summary failed completeness checks", call. = FALSE)
}

components_tsv <- file.path(output_dir, "svclone_native_components_long.tsv")
components_rds <- file.path(output_dir, "svclone_native_components_long.rds")
events_tsv <- file.path(output_dir, "svclone_native_events_long.tsv")
events_rds <- file.path(output_dir, "svclone_native_events_long.rds")
pairwise_tsv <- file.path(output_dir, "svclone_fair_assisted_shared.tsv")
pairwise_rds <- file.path(output_dir, "svclone_fair_assisted_shared.rds")
condition_tsv <- file.path(output_dir, "svclone_arm_condition_summary.tsv")

write.table(components, components_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(components, components_rds, compress = "xz")
write.table(events, events_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(events, events_rds, compress = "xz")
write.table(pairwise, pairwise_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(pairwise, pairwise_rds, compress = "xz")
write.table(condition_summary, condition_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat(sprintf("Assisted: %d native components / %d events across %d conditions\n",
            nrow(assisted), nrow(assisted_events), length(unique(assisted$task_id))))
cat(sprintf("Fair: %d native components / %d events across %d conditions\n",
            nrow(fair), nrow(fair_events), length(unique(fair$task_id))))
cat(sprintf("Fair-assisted shared: %d SVs across %d conditions\n",
            nrow(pairwise), length(unique(pairwise$task_id))))
cat("SVclone arm extraction gate: PASS\n")
