#!/usr/bin/env Rscript

# Join SVCFit to the two SVclone arms without requiring identical breakpoint
# coordinates. Matching is one-to-one within condition and event class, with a
# configurable per-breakpoint tolerance. Partial outputs are deliberately
# prefixed so an in-flight audit cannot be mistaken for the final analysis.

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0L) {
  cat(paste0(
    "Usage: 09_join_three_arms.R --analysis-root DIR ",
    "[--tolerance BP] [--output-prefix TEXT] [--summary-only|--rds-only] [--allow-partial]\n"
  ), file = if (status == 0L) stdout() else stderr())
  quit(status = status)
}

allow_partial <- "--allow-partial" %in% args
summary_only <- "--summary-only" %in% args
rds_only <- "--rds-only" %in% args
args <- args[!args %in% c("--allow-partial", "--summary-only", "--rds-only")]
if (length(args) == 1L && args %in% c("-h", "--help")) usage()
if (length(args) %% 2L) usage(2L)
opt <- list()
for (i in seq(1L, length(args), by = 2L)) {
  key <- sub("^--", "", args[[i]])
  if (identical(key, args[[i]]) || !nzchar(key)) usage(2L)
  opt[[key]] <- args[[i + 1L]]
}
if (is.null(opt[["analysis-root"]])) usage(2L)

root <- normalizePath(opt[["analysis-root"]], mustWork = TRUE)
tolerance <- if (is.null(opt$tolerance)) 100L else suppressWarnings(as.integer(opt$tolerance))
if (!is.finite(tolerance) || tolerance < 0L) stop("--tolerance must be a non-negative integer", call. = FALSE)
if (!requireNamespace("igraph", quietly = TRUE)) stop("The igraph package is required", call. = FALSE)
output_prefix <- if (is.null(opt[["output-prefix"]])) "" else opt[["output-prefix"]]
if (grepl("[/\\\\]", output_prefix) || grepl("^\\.", output_prefix)) {
  stop("--output-prefix must be a filename prefix, not a path", call. = FALSE)
}
if (allow_partial && nzchar(output_prefix)) {
  stop("--output-prefix and --allow-partial cannot be combined", call. = FALSE)
}
if (summary_only && rds_only) stop("--summary-only and --rds-only cannot be combined", call. = FALSE)

manifest_path <- file.path(root, "three_arm_condition_manifest.tsv")
svclone_path <- file.path(root, "svclone_native_events_long.rds")
if (!file.exists(manifest_path) || !file.exists(svclone_path)) {
  stop("Missing manifest or SVclone event extract under analysis root", call. = FALSE)
}
manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
svclone <- readRDS(svclone_path)

svcfit_cache <- file.path(root, "svcfit_events_long.rds")
if (!allow_partial && file.exists(svcfit_cache)) {
  svcfit <- readRDS(svcfit_cache)
} else {
  event_paths <- list.files(file.path(root, "svcfit_events"), pattern = "\\.tsv$",
                            recursive = TRUE, full.names = TRUE)
  event_paths <- event_paths[!grepl("\\.(components|status)\\.tsv$", event_paths)]
  if (!length(event_paths)) stop("No completed SVCFit event tables found", call. = FALSE)
  svcfit_list <- lapply(event_paths, read.delim, stringsAsFactors = FALSE, check.names = FALSE)
  svcfit <- do.call(rbind, svcfit_list)
}
rownames(svcfit) <- NULL

required_svcfit <- c(
  "task_id", "mutation_id", "svcfit_svcf", "true_svcf", "classification",
  "score_ready", "replicate", "experiment", "condition", "truth_purity_pct", "mixture_pct"
)
required_svclone <- c(
  "task_id", "mutation_id", "arm", "chr1", "pos1", "chr2", "pos2",
  "classification", "svcf_native"
)
if (length(setdiff(required_svcfit, names(svcfit)))) stop("Invalid SVCFit event schema", call. = FALSE)
if (length(setdiff(required_svclone, names(svclone)))) stop("Invalid SVclone event schema", call. = FALSE)
if (anyDuplicated(svcfit[c("task_id", "mutation_id")])) stop("Duplicate SVCFit event keys", call. = FALSE)
if (anyDuplicated(svclone[c("task_id", "mutation_id", "arm")])) stop("Duplicate SVclone event keys", call. = FALSE)

completed <- sort(unique(as.integer(svcfit$task_id)))
if (any(!completed %in% manifest$task_id)) stop("SVCFit output contains unknown task IDs", call. = FALSE)
if (!allow_partial && (length(completed) != nrow(manifest) || !identical(completed, sort(manifest$task_id)))) {
  stop(sprintf("Only %d/%d SVCFit conditions are complete; use --allow-partial for an audit",
               length(completed), nrow(manifest)), call. = FALSE)
}
prefix <- if (allow_partial) "partial_" else output_prefix
if (!allow_partial && !file.exists(svcfit_cache)) {
  saveRDS(svcfit, svcfit_cache, compress = "xz")
}

parse_event_id <- function(x) {
  pieces <- strsplit(as.character(x), ":", fixed = TRUE)
  if (any(lengths(pieces) != 4L)) stop("Malformed SVCFit mutation_id", call. = FALSE)
  matrix(unlist(pieces, use.names = FALSE), ncol = 4L, byrow = TRUE,
         dimnames = list(NULL, c("chr1", "pos1", "chr2", "pos2")))
}
coords <- parse_event_id(svcfit$mutation_id)
svcfit$chr1 <- coords[, "chr1"]
svcfit$pos1 <- suppressWarnings(as.integer(coords[, "pos1"]))
svcfit$chr2 <- coords[, "chr2"]
svcfit$pos2 <- suppressWarnings(as.integer(coords[, "pos2"]))
if (anyNA(svcfit[c("chr1", "pos1", "chr2", "pos2")])) stop("Invalid SVCFit coordinates", call. = FALSE)

canonical_class <- function(x) {
  x <- toupper(as.character(x))
  x[x == "BND"] <- "INTRX"
  x
}
svcfit$match_class <- canonical_class(svcfit$classification)
svclone$match_class <- canonical_class(svclone$classification)
svcfit_by_task <- split(svcfit, as.character(svcfit$task_id))
svclone_by_task_arm <- split(
  svclone,
  paste(svclone$task_id, svclone$arm, sep = "\r")
)

# Return a maximum-cardinality, minimum-distance bipartite matching. Candidate
# generation accepts reversed endpoint order. A cardinality-dominating edge
# constant ensures that one additional match is always worth more than every
# possible distance improvement, then total breakpoint distance breaks ties.
match_one <- function(svc, sv, tol) {
  empty <- data.frame(svc_row = integer(), sv_row = integer(), delta1 = integer(),
                      delta2 = integer(), total_delta = integer(), orientation = character(),
                      svc_candidates = integer(), sv_candidates = integer())
  if (!nrow(svc) || !nrow(sv)) return(empty)
  candidates <- merge(
    data.frame(svc_row = seq_len(nrow(svc)), match_class = svc$match_class),
    data.frame(sv_row = seq_len(nrow(sv)), match_class = sv$match_class),
    by = "match_class", sort = FALSE
  )
  if (!nrow(candidates)) return(empty)
  i <- candidates$svc_row
  j <- candidates$sv_row
  direct <- svc$chr1[i] == sv$chr1[j] & svc$chr2[i] == sv$chr2[j]
  reverse <- svc$chr1[i] == sv$chr2[j] & svc$chr2[i] == sv$chr1[j]
  d1_direct <- abs(svc$pos1[i] - sv$pos1[j])
  d2_direct <- abs(svc$pos2[i] - sv$pos2[j])
  d1_reverse <- abs(svc$pos1[i] - sv$pos2[j])
  d2_reverse <- abs(svc$pos2[i] - sv$pos1[j])
  use_reverse <- reverse & (!direct | d1_reverse + d2_reverse < d1_direct + d2_direct)
  candidates$delta1 <- ifelse(use_reverse, d1_reverse, d1_direct)
  candidates$delta2 <- ifelse(use_reverse, d2_reverse, d2_direct)
  candidates$orientation <- ifelse(use_reverse, "reversed", "direct")
  candidates <- candidates[(direct | reverse) & candidates$delta1 <= tol & candidates$delta2 <= tol, ]
  if (!nrow(candidates)) return(empty)
  candidates$total_delta <- candidates$delta1 + candidates$delta2
  svc_n <- table(candidates$svc_row)
  sv_n <- table(candidates$sv_row)
  candidates$svc_candidates <- as.integer(svc_n[as.character(candidates$svc_row)])
  candidates$sv_candidates <- as.integer(sv_n[as.character(candidates$sv_row)])
  candidates <- candidates[order(candidates$total_delta,
                                   pmax(candidates$delta1, candidates$delta2),
                                   candidates$svc_row, candidates$sv_row), ]
  vertex_count <- nrow(svc) + nrow(sv)
  graph <- igraph::make_empty_graph(vertex_count, directed = FALSE)
  graph <- igraph::add_edges(
    graph,
    as.vector(rbind(candidates$svc_row, nrow(svc) + candidates$sv_row))
  )
  max_pairs <- min(nrow(svc), nrow(sv))
  cardinality_constant <- (2L * tol + 1L) * (max_pairs + 1L)
  edge_weight <- cardinality_constant - candidates$total_delta -
    pmax(candidates$delta1, candidates$delta2) / (2L * tol + 1L) -
    seq_len(nrow(candidates)) / ((nrow(candidates) + 1) * 1e6)
  matching <- igraph::max_bipartite_match(
    graph,
    types = c(rep(FALSE, nrow(svc)), rep(TRUE, nrow(sv))),
    weights = edge_weight
  )$matching
  selected_svc <- which(!is.na(matching[seq_len(nrow(svc))]))
  selected_sv <- as.integer(matching[selected_svc] - nrow(svc))
  selected_key <- paste(selected_svc, selected_sv, sep = "\r")
  candidate_key <- paste(candidates$svc_row, candidates$sv_row, sep = "\r")
  selected_rows <- match(selected_key, candidate_key)
  if (anyNA(selected_rows)) stop("Bipartite matching returned a non-candidate edge", call. = FALSE)
  candidates[selected_rows, c("svc_row", "sv_row", "delta1", "delta2", "total_delta",
                               "orientation", "svc_candidates", "sv_candidates")]
}

make_join <- function(task_id, arm) {
  svc <- svcfit_by_task[[as.character(task_id)]]
  sv <- svclone_by_task_arm[[paste(task_id, arm, sep = "\r")]]
  if (is.null(sv)) sv <- svclone[FALSE, , drop = FALSE]
  selected <- match_one(svc, sv, tolerance)
  if (!nrow(selected)) return(NULL)
  a <- svc[selected$svc_row, , drop = FALSE]
  b <- sv[selected$sv_row, , drop = FALSE]
  data.frame(
    task_id = task_id, replicate = a$replicate, experiment = a$experiment,
    condition = a$condition, truth_purity_pct = a$truth_purity_pct,
    mixture_pct = a$mixture_pct, arm = arm,
    svcfit_mutation_id = a$mutation_id, svclone_mutation_id = b$mutation_id,
    classification = a$classification, match_class = a$match_class,
    svcfit_chr1 = a$chr1, svcfit_pos1 = a$pos1,
    svcfit_chr2 = a$chr2, svcfit_pos2 = a$pos2,
    svclone_chr1 = b$chr1, svclone_pos1 = b$pos1,
    svclone_chr2 = b$chr2, svclone_pos2 = b$pos2,
    breakpoint_delta1 = selected$delta1, breakpoint_delta2 = selected$delta2,
    total_breakpoint_delta = selected$total_delta,
    match_orientation = selected$orientation,
    svcfit_candidate_count = selected$svc_candidates,
    svclone_candidate_count = selected$sv_candidates,
    match_type = ifelse(selected$total_delta == 0L, "exact", "tolerant"),
    true_svcf = a$true_svcf, svcfit_svcf = a$svcfit_svcf,
    svclone_svcf = b$svcf_native,
    svcfit_signed_error = a$true_svcf - a$svcfit_svcf,
    svcfit_abs_error = abs(a$true_svcf - a$svcfit_svcf),
    svclone_signed_error = a$true_svcf - b$svcf_native,
    svclone_abs_error = abs(a$true_svcf - b$svcf_native),
    score_ready = as.logical(a$score_ready) & is.finite(b$svcf_native),
    stringsAsFactors = FALSE
  )
}

joined <- list()
n <- 0L
for (task_id in completed) {
  for (arm in c("assisted", "fair")) {
    x <- make_join(task_id, arm)
    if (!is.null(x)) {
      n <- n + 1L
      joined[[n]] <- x
    }
  }
}
if (!length(joined)) stop("No SVCFit/SVclone matches found", call. = FALSE)
joined <- do.call(rbind, joined)
rownames(joined) <- NULL
if (anyDuplicated(joined[c("task_id", "arm", "svcfit_mutation_id")]) ||
    anyDuplicated(joined[c("task_id", "arm", "svclone_mutation_id")])) {
  stop("One-to-one matching invariant failed", call. = FALSE)
}

count_task_arm <- function(x, value_name) {
  tab <- as.data.frame(table(task_id = x$task_id, arm = x$arm), stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0L, ]
  names(tab)[names(tab) == "Freq"] <- value_name
  tab$task_id <- as.integer(as.character(tab$task_id))
  tab
}
svc_counts <- data.frame(task_id = as.integer(names(table(svcfit$task_id))),
                         svcfit_event_count = as.integer(table(svcfit$task_id)))
sv_counts <- count_task_arm(svclone[svclone$task_id %in% completed, ], "svclone_event_count")
match_counts <- count_task_arm(joined, "matched_event_count")
condition_summary <- merge(
  expand.grid(task_id = completed, arm = c("assisted", "fair"), stringsAsFactors = FALSE),
  svc_counts, by = "task_id", all.x = TRUE, sort = FALSE
)
condition_summary <- merge(condition_summary, sv_counts, by = c("task_id", "arm"), all.x = TRUE, sort = FALSE)
condition_summary <- merge(condition_summary, match_counts, by = c("task_id", "arm"), all.x = TRUE, sort = FALSE)
for (nm in c("svclone_event_count", "matched_event_count")) condition_summary[[nm]][is.na(condition_summary[[nm]])] <- 0L
condition_summary$svcfit_match_fraction <- condition_summary$matched_event_count / condition_summary$svcfit_event_count
condition_summary$svclone_match_fraction <- ifelse(condition_summary$svclone_event_count > 0L,
                                                    condition_summary$matched_event_count / condition_summary$svclone_event_count,
                                                    NA_real_)
condition_summary <- merge(condition_summary,
                           manifest[c("task_id", "replicate", "experiment", "condition",
                                      "truth_purity_pct", "mixture_pct", "fair_outcome")],
                           by = "task_id", all.x = TRUE, sort = FALSE)
condition_summary <- condition_summary[order(condition_summary$task_id, condition_summary$arm), ]

score <- joined[joined$score_ready, ]
accuracy <- do.call(rbind, lapply(split(score, score$arm), function(x) data.frame(
  arm = x$arm[[1L]], conditions = length(unique(x$task_id)), events = nrow(x),
  svcfit_mae = mean(x$svcfit_abs_error), svcfit_bias = mean(x$svcfit_signed_error),
  svclone_mae = mean(x$svclone_abs_error), svclone_bias = mean(x$svclone_signed_error),
  exact_matches = sum(x$match_type == "exact"), tolerant_matches = sum(x$match_type == "tolerant"),
  ambiguous_matches = sum(x$svcfit_candidate_count > 1L | x$svclone_candidate_count > 1L)
)))
rownames(accuracy) <- NULL

assisted <- joined[joined$arm == "assisted", ]
fair <- joined[joined$arm == "fair", ]
three_way <- merge(
  assisted, fair,
  by = c("task_id", "svcfit_mutation_id"), suffixes = c("_assisted", "_fair"),
  all = FALSE, sort = FALSE
)

if (!summary_only) {
  if (!rds_only) {
    write.table(joined, file.path(root, paste0(prefix, "three_arm_pairwise_matches.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  saveRDS(joined, file.path(root, paste0(prefix, "three_arm_pairwise_matches.rds")), compress = "xz")
  if (!rds_only) {
    write.table(three_way, file.path(root, paste0(prefix, "three_arm_shared_matches.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  saveRDS(three_way, file.path(root, paste0(prefix, "three_arm_shared_matches.rds")), compress = "xz")
}
write.table(condition_summary, file.path(root, paste0(prefix, "three_arm_match_by_condition.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(accuracy, file.path(root, paste0(prefix, "three_arm_accuracy_summary.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat(sprintf("SVCFit conditions: %d/%d%s\n", length(completed), nrow(manifest),
            if (allow_partial) " (partial audit)" else ""))
cat(sprintf("Pairwise matches: %d; three-way shared: %d; tolerance: %d bp\n",
            nrow(joined), nrow(three_way), tolerance))
print(accuracy, row.names = FALSE)
cat("Three-arm join gate: PASS\n")
