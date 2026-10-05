#!/usr/bin/env Rscript

# Summarize the validated three-arm matches without treating event rows as
# independent replicates. Point estimates retain event weighting; percentile
# confidence intervals resample whole simulation replicates.

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0L) {
  cat(paste0(
    "Usage: 10_summarize_three_arms.R --analysis-root DIR ",
    "[--output-dir DIR] [--bootstrap N] [--seed N]\n"
  ), file = if (status == 0L) stdout() else stderr())
  quit(status = status)
}
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
output_dir <- if (is.null(opt[["output-dir"]])) {
  file.path(root, "three_arm_summary")
} else {
  normalizePath(opt[["output-dir"]], mustWork = FALSE)
}
bootstrap_n <- if (is.null(opt$bootstrap)) 1000L else suppressWarnings(as.integer(opt$bootstrap))
seed <- if (is.null(opt$seed)) 20260920L else suppressWarnings(as.integer(opt$seed))
if (!is.finite(bootstrap_n) || bootstrap_n < 100L) stop("--bootstrap must be at least 100", call. = FALSE)
if (!is.finite(seed)) stop("--seed must be an integer", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory", call. = FALSE)

paths <- c(
  matches = file.path(root, "three_arm_pairwise_matches.rds"),
  conditions = file.path(root, "three_arm_match_by_condition.tsv"),
  manifest = file.path(root, "three_arm_condition_manifest.tsv")
)
if (any(!file.exists(paths))) stop("Missing final three-arm inputs", call. = FALSE)
matches <- readRDS(paths[["matches"]])
conditions <- read.delim(paths[["conditions"]], stringsAsFactors = FALSE, check.names = FALSE)
manifest <- read.delim(paths[["manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

required_matches <- c(
  "task_id", "replicate", "experiment", "condition", "truth_purity_pct",
  "mixture_pct", "arm", "svcfit_mutation_id", "match_class", "true_svcf",
  "svcfit_svcf", "svclone_svcf", "score_ready"
)
if (length(setdiff(required_matches, names(matches)))) stop("Invalid final match schema", call. = FALSE)
if (nrow(matches) != 275326L || anyDuplicated(matches[c("task_id", "arm", "svcfit_mutation_id")])) {
  stop("Final pairwise match gate failed", call. = FALSE)
}
if (nrow(manifest) != 2250L || sum(manifest$fair_outcome == "ready") != 1727L ||
    sum(manifest$fair_outcome == "failed_no_purity") != 523L) {
  stop("Manifest availability gate failed", call. = FALSE)
}

metadata <- c(
  "task_id", "replicate", "experiment", "condition", "truth_purity_pct",
  "mixture_pct", "svcfit_mutation_id", "match_class", "true_svcf"
)
as_method <- function(x, scope, method, estimate_column) {
  out <- x[metadata]
  out$scope <- scope
  out$method <- method
  out$estimate <- as.numeric(x[[estimate_column]])
  out$signed_error <- out$true_svcf - out$estimate
  out$abs_error <- abs(out$signed_error)
  out$sq_error <- out$signed_error^2
  out[is.finite(out$true_svcf) & is.finite(out$estimate), ]
}

assisted <- matches[matches$arm == "assisted" & matches$score_ready, ]
fair <- matches[matches$arm == "fair" & matches$score_ready, ]
assisted_key <- paste(assisted$task_id, assisted$svcfit_mutation_id, sep = "\r")
fair_key <- paste(fair$task_id, fair$svcfit_mutation_id, sep = "\r")
three_key <- intersect(assisted_key, fair_key)
assisted_three <- assisted[match(three_key, assisted_key), ]
fair_three <- fair[match(three_key, fair_key), ]
if (anyNA(assisted_three$task_id) || anyNA(fair_three$task_id) ||
    any(abs(assisted_three$true_svcf - fair_three$true_svcf) > 1e-12) ||
    any(abs(assisted_three$svcfit_svcf - fair_three$svcfit_svcf) > 1e-12)) {
  stop("Three-way shared-set alignment failed", call. = FALSE)
}

evaluation <- rbind(
  as_method(assisted, "pairwise_assisted", "SVCFit", "svcfit_svcf"),
  as_method(assisted, "pairwise_assisted", "SVclone assisted", "svclone_svcf"),
  as_method(fair, "pairwise_fair", "SVCFit", "svcfit_svcf"),
  as_method(fair, "pairwise_fair", "SVclone fair", "svclone_svcf"),
  as_method(assisted_three, "three_way", "SVCFit", "svcfit_svcf"),
  as_method(assisted_three, "three_way", "SVclone assisted", "svclone_svcf"),
  as_method(fair_three, "three_way", "SVclone fair", "svclone_svcf")
)
rownames(evaluation) <- NULL

set.seed(seed)
cluster_ci <- function(x) {
  clusters <- split(seq_len(nrow(x)), x$replicate)
  k <- length(clusters)
  if (k < 2L) return(c(mae_low = NA_real_, mae_high = NA_real_,
                        bias_low = NA_real_, bias_high = NA_real_))
  cluster_n <- vapply(clusters, length, integer(1))
  cluster_abs <- vapply(clusters, function(i) sum(x$abs_error[i]), numeric(1))
  cluster_signed <- vapply(clusters, function(i) sum(x$signed_error[i]), numeric(1))
  sampled <- matrix(sample.int(k, k * bootstrap_n, replace = TRUE), nrow = k)
  denom <- colSums(matrix(cluster_n[sampled], nrow = k))
  mae <- colSums(matrix(cluster_abs[sampled], nrow = k)) / denom
  bias <- colSums(matrix(cluster_signed[sampled], nrow = k)) / denom
  c(
    mae_low = unname(quantile(mae, 0.025, names = FALSE)),
    mae_high = unname(quantile(mae, 0.975, names = FALSE)),
    bias_low = unname(quantile(bias, 0.025, names = FALSE)),
    bias_high = unname(quantile(bias, 0.975, names = FALSE))
  )
}

summarize_accuracy <- function(x, scope, method, stratum, value) {
  ci <- cluster_ci(x)
  data.frame(
    scope = scope, method = method, stratum = stratum,
    stratum_value = as.character(value), events = nrow(x),
    conditions = length(unique(x$task_id)), replicates = length(unique(x$replicate)),
    mae = mean(x$abs_error), mae_ci_low = ci[["mae_low"]], mae_ci_high = ci[["mae_high"]],
    rmse = sqrt(mean(x$sq_error)), median_abs_error = median(x$abs_error),
    bias = mean(x$signed_error), bias_ci_low = ci[["bias_low"]],
    bias_ci_high = ci[["bias_high"]], stringsAsFactors = FALSE
  )
}

strata <- list(
  overall = function(x) rep("all", nrow(x)),
  truth_purity_pct = function(x) x$truth_purity_pct,
  mixture_pct = function(x) x$mixture_pct,
  experiment = function(x) x$experiment,
  condition = function(x) x$condition,
  sv_class = function(x) x$match_class
)
accuracy_rows <- list()
row_index <- 0L
evaluation_groups <- split(seq_len(nrow(evaluation)),
                           paste(evaluation$scope, evaluation$method, sep = "\r"))
for (idx in evaluation_groups) {
  x <- evaluation[idx, , drop = FALSE]
  for (stratum in names(strata)) {
    values <- as.character(strata[[stratum]](x))
    for (value in sort(unique(values))) {
      row_index <- row_index + 1L
      accuracy_rows[[row_index]] <- summarize_accuracy(
        x[values == value, , drop = FALSE], x$scope[[1L]], x$method[[1L]], stratum, value
      )
    }
  }
}
accuracy <- do.call(rbind, accuracy_rows)
rownames(accuracy) <- NULL
accuracy <- accuracy[order(accuracy$scope, accuracy$stratum,
                           accuracy$stratum_value, accuracy$method), ]

condition_strata <- list(
  overall = function(x) rep("all", nrow(x)),
  truth_purity_pct = function(x) x$truth_purity_pct,
  mixture_pct = function(x) x$mixture_pct,
  experiment = function(x) x$experiment,
  condition = function(x) x$condition
)
coverage_rows <- list()
row_index <- 0L
for (arm in c("assisted", "fair")) {
  arm_data <- conditions[conditions$arm == arm, ]
  for (stratum in names(condition_strata)) {
    values <- as.character(condition_strata[[stratum]](arm_data))
    for (value in sort(unique(values))) {
      z <- arm_data[values == value, , drop = FALSE]
      available <- z$svclone_event_count > 0L
      row_index <- row_index + 1L
      coverage_rows[[row_index]] <- data.frame(
        arm = arm, stratum = stratum, stratum_value = value,
        conditions_total = nrow(z), conditions_available = sum(available),
        svcfit_events = sum(z$svcfit_event_count),
        svclone_events = sum(z$svclone_event_count), matched_events = sum(z$matched_event_count),
        svcfit_event_weighted_match_fraction = sum(z$matched_event_count) / sum(z$svcfit_event_count),
        svclone_event_weighted_match_fraction = if (sum(z$svclone_event_count) > 0L)
          sum(z$matched_event_count) / sum(z$svclone_event_count) else NA_real_,
        median_svcfit_match_fraction = if (any(available)) median(z$svcfit_match_fraction[available]) else NA_real_,
        median_svclone_match_fraction = if (any(available)) median(z$svclone_match_fraction[available]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}
coverage <- do.call(rbind, coverage_rows)
rownames(coverage) <- NULL

availability_rows <- list()
row_index <- 0L
for (stratum in names(condition_strata)) {
  values <- as.character(condition_strata[[stratum]](manifest))
  for (value in sort(unique(values))) {
    z <- manifest[values == value, , drop = FALSE]
    row_index <- row_index + 1L
    availability_rows[[row_index]] <- data.frame(
      stratum = stratum, stratum_value = value, conditions = nrow(z),
      fair_ready = sum(z$fair_outcome == "ready"),
      fair_failed_no_purity = sum(z$fair_outcome == "failed_no_purity"),
      fair_ready_fraction = mean(z$fair_outcome == "ready"), stringsAsFactors = FALSE
    )
  }
}
availability <- do.call(rbind, availability_rows)
rownames(availability) <- NULL

read_sensitivity <- function(tolerance) {
  prefix <- if (tolerance == 100L) "" else sprintf("sensitivity_tol%03d_", tolerance)
  accuracy_path <- file.path(root, paste0(prefix, "three_arm_accuracy_summary.tsv"))
  condition_path <- file.path(root, paste0(prefix, "three_arm_match_by_condition.tsv"))
  if (!file.exists(accuracy_path) || !file.exists(condition_path)) return(NULL)
  a <- read.delim(accuracy_path, stringsAsFactors = FALSE, check.names = FALSE)
  cnd <- read.delim(condition_path, stringsAsFactors = FALSE, check.names = FALSE)
  totals <- do.call(rbind, lapply(split(cnd, cnd$arm), function(z) data.frame(
    arm = z$arm[[1L]], svcfit_events = sum(z$svcfit_event_count),
    svclone_events = sum(z$svclone_event_count), matched_events = sum(z$matched_event_count),
    svcfit_match_fraction = sum(z$matched_event_count) / sum(z$svcfit_event_count),
    svclone_match_fraction = sum(z$matched_event_count) / sum(z$svclone_event_count)
  )))
  out <- merge(a, totals, by = "arm", sort = FALSE)
  out$tolerance_bp <- tolerance
  out
}
sensitivity <- do.call(rbind, Filter(Negate(is.null), lapply(c(0L, 25L, 50L, 100L), read_sensitivity)))
rownames(sensitivity) <- NULL
sensitivity <- sensitivity[order(sensitivity$arm, sensitivity$tolerance_bp), ]

write.table(accuracy, file.path(output_dir, "three_arm_accuracy_stratified.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(coverage, file.path(output_dir, "three_arm_coverage_stratified.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(availability, file.path(output_dir, "three_arm_fair_availability_stratified.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(sensitivity, file.path(output_dir, "three_arm_tolerance_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(evaluation, file.path(output_dir, "three_arm_evaluation_long.rds"), compress = "xz")

method_colors <- c("SVCFit" = "#31688e", "SVclone assisted" = "#35b779", "SVclone fair" = "#d95f02")
overall <- accuracy[accuracy$scope == "three_way" & accuracy$stratum == "overall", ]
overall <- overall[match(names(method_colors), overall$method), ]
pdf(file.path(output_dir, "three_arm_overall_mae.pdf"), width = 7, height = 5)
bp <- barplot(overall$mae, names.arg = overall$method, col = method_colors[overall$method],
              ylab = "Mean absolute SVCF error", ylim = c(0, max(overall$mae_ci_high) * 1.15),
              main = "Three-way shared events")
arrows(bp, overall$mae_ci_low, bp, overall$mae_ci_high, angle = 90, code = 3, length = 0.05)
mtext(sprintf("Replicate bootstrap: %d resamples", bootstrap_n), side = 1, line = 4, cex = 0.8)
dev.off()

purity <- accuracy[accuracy$scope == "three_way" & accuracy$stratum == "truth_purity_pct", ]
purity$purity <- as.numeric(purity$stratum_value)
pdf(file.path(output_dir, "three_arm_mae_by_purity.pdf"), width = 7, height = 5)
plot(range(purity$purity), range(c(purity$mae_ci_low, purity$mae_ci_high)), type = "n",
     xlab = "Truth purity (%)", ylab = "Mean absolute SVCF error",
     main = "Three-way accuracy by truth purity")
for (method in names(method_colors)) {
  z <- purity[purity$method == method, ]
  z <- z[order(z$purity), ]
  lines(z$purity, z$mae, type = "b", pch = 19, col = method_colors[[method]])
  arrows(z$purity, z$mae_ci_low, z$purity, z$mae_ci_high,
         angle = 90, code = 3, length = 0.03, col = method_colors[[method]])
}
legend("topleft", names(method_colors), col = method_colors, lty = 1, pch = 19, bty = "n")
dev.off()

pdf(file.path(output_dir, "three_arm_tolerance_sensitivity.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(NA, xlim = range(sensitivity$tolerance_bp), ylim = range(sensitivity$events),
     xlab = "Breakpoint tolerance (bp)", ylab = "Matched events", main = "Event retention")
for (arm in c("assisted", "fair")) {
  z <- sensitivity[sensitivity$arm == arm, ]
  lines(z$tolerance_bp, z$events, type = "b", pch = 19,
        col = if (arm == "assisted") method_colors[["SVclone assisted"]] else method_colors[["SVclone fair"]])
}
legend("bottomright", c("Assisted", "Fair"),
       col = c(method_colors[["SVclone assisted"]], method_colors[["SVclone fair"]]),
       lty = 1, pch = 19, bty = "n")
plot(NA, xlim = range(sensitivity$tolerance_bp),
     ylim = range(c(sensitivity$svcfit_mae, sensitivity$svclone_mae)),
     xlab = "Breakpoint tolerance (bp)", ylab = "Mean absolute SVCF error",
     main = "Pairwise-set MAE")
for (arm in c("assisted", "fair")) {
  z <- sensitivity[sensitivity$arm == arm, ]
  arm_color <- if (arm == "assisted") method_colors[["SVclone assisted"]] else method_colors[["SVclone fair"]]
  lines(z$tolerance_bp, z$svclone_mae, type = "b", pch = 19, col = arm_color)
  lines(z$tolerance_bp, z$svcfit_mae, type = "b", pch = 1, lty = 2, col = arm_color)
}
legend("topright", c("SVclone assisted", "SVCFit on assisted set", "SVclone fair", "SVCFit on fair set"),
       col = c(method_colors[["SVclone assisted"]], method_colors[["SVclone assisted"]],
               method_colors[["SVclone fair"]], method_colors[["SVclone fair"]]),
       lty = c(1, 2, 1, 2), pch = c(19, 1, 19, 1), cex = 0.75, bty = "n")
dev.off()

cat(sprintf("Evaluation rows: %d; three-way events: %d\n", nrow(evaluation), length(three_key)))
cat(sprintf("Accuracy strata: %d rows; bootstrap replicates: %d\n", nrow(accuracy), bootstrap_n))
cat(sprintf("Sensitivity tolerances: %s\n", paste(sort(unique(sensitivity$tolerance_bp)), collapse = ", ")))
cat("Three-arm summary gate: PASS\n")
