#!/usr/bin/env Rscript

# Build the prostate-mixture comparison on the identical detected-SV set:
#   1. SVCFit
#   2. SVclone truth-CCF-assisted (known CCF selects multiplicity)
#   3. SVclone non-assisted (native CCube CCF)
#
# This script consumes completed immutable runs. It does not rerun callers or
# overwrite their outputs. All machine-specific paths are required arguments.

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: 08_prostate_three_arm.R \\\n",
    "  --svcfit-chrx FILE --svcfit-45 FILE \\\n",
    "  --assisted-root DIR --fair-root DIR --truth-dir DIR \\\n",
    "  --output-dir DIR --svclone-rlib DIR [--bootstrap 2000] [--seed 20260920]\n"
  ))
  quit(status = status)
}

parse_args <- function(x) {
  if (any(x %in% c("-h", "--help"))) usage(0L)
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!grepl("^--", x[[i]]) || i == length(x)) usage(2L)
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x
opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("svcfit-chrx", "svcfit-45", "assisted-root", "fair-root",
              "truth-dir", "output-dir", "svclone-rlib")
missing <- required[!vapply(required, function(n) nzchar(opt[[n]] %||% ""), logical(1))]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "), call. = FALSE)

bootstrap_n <- as.integer(opt[["bootstrap"]] %||% "2000")
seed <- as.integer(opt[["seed"]] %||% "20260920")
if (!is.finite(bootstrap_n) || bootstrap_n < 100L) stop("--bootstrap must be >= 100")

inputs <- unlist(opt[required[required != "output-dir"]], use.names = FALSE)
absent <- inputs[!file.exists(inputs)]
if (length(absent)) stop("Missing inputs:\n", paste(absent, collapse = "\n"), call. = FALSE)
dir.create(opt[["output-dir"]], recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(opt[["output-dir"]])

.libPaths(c(normalizePath(opt[["svclone-rlib"]]), .libPaths()))
suppressPackageStartupMessages({
  library(ccube)
  library(dplyr)
  library(stringr)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
helper_dir <- file.path(script_dir, "helper")
source(file.path(helper_dir, "util.R"))

conditions <- c("3m19", "3m28", "3m37", "3m46", "3m55",
                "3m64", "3m73", "3m82", "3m91", "4m", "5m")
reps <- paste0("rep", seq_len(30L))

# The original helper rebuilt this mapping, including an O(n_bM * n_gM)
# overlap search, for every sample. Truth is invariant across all 330 samples,
# so build the exact coordinate-to-source map once.
read_truth_ids <- function(path) {
  x <- read.table(path, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  names(x) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                "INFO", "FORMAT", "normal", "tumor")
  x <- x[x$FILTER == "PASS", ]
  is_bnd <- grepl("BND", x$ID)
  x$pos2 <- ifelse(is_bnd, str_extract(x$ALT, "(?<=:)\\d+"),
                   sub(".*END=(\\d+);.*", "\\1", x$INFO))
  x$chr2 <- ifelse(is_bnd, str_extract(x$ALT, "(?<=\\[|\\])[^:]+"), x$CHROM)
  x$POS <- as.numeric(x$POS); x$pos2 <- as.numeric(x$pos2)
  x <- x[as.character(x$CHROM) %in% c(as.character(1:22), "X") &
           as.character(x$chr2) %in% c(as.character(1:22), "X"), ]
  x$id <- paste(x$CHROM, x$POS, x$chr2, x$pos2, sep = ":")
  x
}

bm_truth <- read_truth_ids(file.path(opt[["truth-dir"]], "bM.vcf"))
gm_truth <- read_truth_ids(file.path(opt[["truth-dir"]], "gM.vcf"))
bm_shared <- logical(nrow(bm_truth)); gm_shared <- logical(nrow(gm_truth))
for (i in seq_len(nrow(bm_truth))) {
  hit <- which(as.character(gm_truth$CHROM) == as.character(bm_truth$CHROM[[i]]) &
               as.character(gm_truth$chr2) == as.character(bm_truth$chr2[[i]]) &
               abs(gm_truth$POS - bm_truth$POS[[i]]) < 15 &
               abs(gm_truth$pos2 - bm_truth$pos2[[i]]) < 15)
  if (length(hit)) { bm_shared[[i]] <- TRUE; gm_shared[hit] <- TRUE }
}
truth_source <- c(setNames(ifelse(bm_shared, "shared", "bM_only"), bm_truth$id),
                  setNames(ifelse(gm_shared, "shared", "gM_only"), gm_truth$id))
# Where the same coordinate occurs in both VCFs, shared must win.
truth_source <- tapply(unname(truth_source), names(truth_source),
                       function(z) if ("shared" %in% z) "shared" else z[[1]])

read_svcfit <- function() {
  a <- read.delim(opt[["svcfit-chrx"]], check.names = FALSE, stringsAsFactors = FALSE)
  b <- read.delim(opt[["svcfit-45"]], check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(a), names(b))) stop("SVCFit tables have different schemas")
  x <- bind_rows(a, b)
  needed <- c("CHROM", "POS", "chr2", "END", "final_svcf", "purity", "replicate", "condition")
  if (length(setdiff(needed, names(x)))) stop("SVCFit table missing required columns")
  x$rep <- paste0("rep", as.integer(x$replicate))
  x$condition <- as.character(x$condition)
  x$mutation_id <- paste(sub("^chr", "", x$CHROM), as.integer(x$POS),
                         sub("^chr", "", x$chr2), as.integer(x$END), sep = ":")
  x$svcfit_ccf <- pmin(as.numeric(x$final_svcf) / as.numeric(x$purity), 1)
  x <- x[is.finite(x$svcfit_ccf) & !is.na(x$mutation_id), ]
  split(x, paste(x$rep, x$condition, sep = "\r"))
}

load_ssm <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("doubleBreakPtsRes", envir = e, inherits = FALSE) ||
      is.null(e$doubleBreakPtsRes$ssm)) stop("Invalid SVclone RData: ", path)
  as.data.frame(e$doubleBreakPtsRes$ssm, stringsAsFactors = FALSE)
}

plain_id <- function(x) sub("(.*):[+-]_(.*):[+-]$", "\\1:\\2", x)

truth_ccf <- function(x, condition) {
  ans <- rep(1, nrow(x))
  if (grepl("^3m", condition)) {
    p <- as.integer(strsplit(sub("^3m", "", condition), "")[[1]]) / 10
    ans[x$sample == "bM_only"] <- p[[1]]
    ans[x$sample == "gM_only"] <- p[[2]]
  } else {
    chr <- suppressWarnings(as.integer(as.character(x$chr1)))
    if (condition == "4m") {
      ans[x$sample == "bM_only" & chr %% 2 == 0] <- 0.6
      ans[x$sample == "bM_only" & chr %% 2 == 1] <- 0.2
      ans[x$sample == "gM_only" & chr %% 2 == 1] <- 0.4
      ans[x$sample == "gM_only" & chr %% 2 == 0] <- NA_real_
    } else {
      ans[x$sample == "bM_only" & chr %% 2 == 0] <- 0.6
      ans[x$sample == "bM_only" & chr %% 2 == 1] <- 0.8
      ans[x$sample == "gM_only" & chr %% 2 == 1] <- 0.2
      ans[x$sample == "gM_only" & chr %% 2 == 0] <- 0.4
    }
  }
  ans
}

mean_finite <- function(x) {
  x <- as.numeric(x); if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
}

collapse_method <- function(x, value, truth = NULL) {
  if (!nrow(x)) return(x)
  keys <- unique(x$mutation_id)
  data.frame(
    mutation_id = keys,
    estimate = vapply(keys, function(k) mean_finite(x[[value]][x$mutation_id == k]), numeric(1)),
    truth = if (is.null(truth)) NA_real_ else
      vapply(keys, function(k) mean_finite(x[[truth]][x$mutation_id == k]), numeric(1)),
    stringsAsFactors = FALSE
  )
}

load_assisted <- function(rep, condition) {
  root <- file.path(opt[["assisted-root"]], rep, "svclone", condition, condition)
  rdata <- file.path(root, "ccube_out", paste0(condition, "_ccube_sv_results.RData"))
  filt <- file.path(root, paste0(condition, "_filtered_svs.tsv"))
  pp_file <- file.path(root, "purity_ploidy.txt")
  for (f in c(rdata, filt, pp_file)) if (!file.exists(f)) stop("Missing assisted input: ", f)
  ssm <- load_ssm(rdata)
  ssm$mutation_id <- plain_id(ssm$mutation_id)
  svin <- read.delim(filt, stringsAsFactors = FALSE)
  svin$mutation_id <- paste(svin$chr1, svin$pos1, svin$chr2, svin$pos2, sep = ":")
  svin$sv <- svin$mutation_id
  # SVclone may represent one inversion with multiple signed components, and
  # its filtered table can repeat that coordinate. The Cartesian rows are
  # intentional here and are collapsed to one strand-free event below.
  x <- inner_join(ssm, svin, by = "mutation_id", relationship = "many-to-many")
  x$sample <- unname(truth_source[x$sv])
  x <- x[!is.na(x$sample), ]
  x$true_ccf <- truth_ccf(x, condition)
  x <- x[is.finite(x$true_ccf), ]
  pp <- read.delim(pp_file, stringsAsFactors = FALSE)
  pur <- as.numeric(pp$purity[[1]]); pl <- as.numeric(pp$ploidy[[1]])
  if (!is.finite(pur) || !is.finite(pl)) stop("Invalid assisted purity/ploidy: ", pp_file)
  af <- 1 - pur / pl
  gain <- x$classification %in% c("DUP", "INTDUP")
  x$norm1_adjusted <- x$norm1; x$norm2_adjusted <- x$norm2
  x$norm1_adjusted[gain] <- x$norm1_adjusted[gain] * af
  x$norm2_adjusted[gain] <- x$norm2_adjusted[gain] * af
  x$norm_cn <- ifelse(as.character(x$chr1) %in% c("X", "Y", "MT"), 1, 2)
  a <- t(vapply(seq_len(nrow(x)), function(i)
    get_true_sc_cn_by_side(x[i, ], pur = pur, side = 1), numeric(4)))
  b <- t(vapply(seq_len(nrow(x)), function(i)
    get_true_sc_cn_by_side(x[i, ], pur = pur, side = 2), numeric(4)))
  x$assisted_ccf <- pmin(rowMeans(cbind(a[, 4], b[, 4])), 1)
  z <- collapse_method(x, "assisted_ccf", "true_ccf")
  names(z)[names(z) == "estimate"] <- "assisted_ccf"
  list(data = z, purity = pur, raw_n = nrow(x))
}

load_fair <- function(rep, condition) {
  root <- file.path(opt[["fair-root"]], rep, "svclone", condition, condition)
  rdata <- file.path(root, "ccube_out", paste0(condition, "_ccube_sv_results.RData"))
  if (!file.exists(rdata)) stop("Missing fair input: ", rdata)
  x <- load_ssm(rdata)
  x$mutation_id <- plain_id(x$mutation_id)
  x$fair_ccf <- pmin(rowMeans(cbind(as.numeric(x$ccube_ccf1),
                                    as.numeric(x$ccube_ccf2))), 1)
  z <- collapse_method(x, "fair_ccf")
  z$truth <- NULL
  names(z)[names(z) == "estimate"] <- "fair_ccf"
  list(data = z, purity = mean(as.numeric(x$purity)), raw_n = nrow(x))
}

svcfit <- read_svcfit()
manifest <- vector("list", length(reps) * length(conditions))
matched <- vector("list", length(manifest))
idx <- 0L
for (rep in reps) for (condition in conditions) {
  idx <- idx + 1L
  key <- paste(rep, condition, sep = "\r")
  sf0 <- svcfit[[key]]
  if (is.null(sf0)) stop("No SVCFit rows for ", rep, "/", condition)
  sf <- collapse_method(sf0, "svcfit_ccf")
  sf$truth <- NULL
  names(sf)[names(sf) == "estimate"] <- "svcfit_ccf"
  assisted <- load_assisted(rep, condition)
  fair <- load_fair(rep, condition)
  if (abs(assisted$purity - fair$purity) > 1e-8)
    stop("Purity differs between SVclone arms for ", rep, "/", condition)
  x <- inner_join(assisted$data, fair$data, by = "mutation_id")
  x <- inner_join(x, sf, by = "mutation_id")
  x <- x[is.finite(x$truth) & is.finite(x$assisted_ccf) &
           is.finite(x$fair_ccf) & is.finite(x$svcfit_ccf), ]
  if (!nrow(x)) stop("No three-way matched events for ", rep, "/", condition)
  x$rep <- rep; x$replicate <- as.integer(sub("^rep", "", rep)); x$condition <- condition
  x$error_svcfit <- x$truth - x$svcfit_ccf
  x$error_assisted <- x$truth - x$assisted_ccf
  x$error_fair <- x$truth - x$fair_ccf
  matched[[idx]] <- x
  manifest[[idx]] <- data.frame(
    rep = rep, replicate = x$replicate[[1]], condition = condition,
    svcfit_events = nrow(sf), assisted_events = nrow(assisted$data),
    fair_events = nrow(fair$data), three_way_events = nrow(x),
    assisted_purity = assisted$purity, fair_purity = fair$purity,
    stringsAsFactors = FALSE)
  cat(sprintf("[%3d/330] %s/%s: %d shared events\n", idx, rep, condition, nrow(x)))
}

manifest <- bind_rows(manifest)
matched <- bind_rows(matched)
if (nrow(manifest) != 330L || length(unique(paste(manifest$rep, manifest$condition))) != 330L)
  stop("Manifest gate failed: expected 330 unique conditions")
if (any(manifest$three_way_events <= 0L)) stop("Manifest gate failed: empty shared set")
if (anyDuplicated(matched[c("rep", "condition", "mutation_id")]))
  stop("Matched-event gate failed: duplicate event keys")

methods <- c("SVCFit", "SVclone truth-CCF-assisted", "SVclone non-assisted")
long <- bind_rows(
  data.frame(matched[c("rep", "replicate", "condition", "mutation_id", "truth")],
             method = methods[[1]], estimate = matched$svcfit_ccf, error = matched$error_svcfit),
  data.frame(matched[c("rep", "replicate", "condition", "mutation_id", "truth")],
             method = methods[[2]], estimate = matched$assisted_ccf, error = matched$error_assisted),
  data.frame(matched[c("rep", "replicate", "condition", "mutation_id", "truth")],
             method = methods[[3]], estimate = matched$fair_ccf, error = matched$error_fair)
)
long$method <- factor(long$method, levels = methods)

metric_one <- function(d) data.frame(
  n_events = nrow(d), mean_error = mean(d$error),
  abs_mean_error = abs(mean(d$error)), mae = mean(abs(d$error)), rmse = sqrt(mean(d$error^2)))
per_rep <- bind_rows(lapply(split(long, interaction(long$rep, long$condition, long$method, drop = TRUE)),
                            metric_one), .id = "group")
parts <- do.call(rbind, strsplit(as.character(per_rep$group), "\\."))
# Conditions and method labels contain no dots; rep is first and method is last.
per_rep$rep <- parts[, 1]
per_rep$condition <- parts[, 2]
per_rep$method <- factor(parts[, 3], levels = methods)
per_rep$replicate <- as.integer(sub("^rep", "", per_rep$rep))
per_rep$group <- NULL

boot_ci <- function(x, B = bootstrap_n) {
  x <- x[is.finite(x)]
  b <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x), lo = unname(quantile(b, .025)), hi = unname(quantile(b, .975)))
}
set.seed(seed)
summary_rows <- list(); j <- 0L
for (condition in conditions) for (method in methods) {
  d <- per_rep[per_rep$condition == condition & per_rep$method == method, ]
  for (metric in c("mean_error", "abs_mean_error", "mae", "rmse")) {
    j <- j + 1L; ci <- boot_ci(d[[metric]])
    summary_rows[[j]] <- data.frame(condition = condition, method = method, metric = metric,
                                    estimate = ci[[1]], ci_lo = ci[[2]], ci_hi = ci[[3]],
                                    n_replicates = nrow(d), stringsAsFactors = FALSE)
  }
}
summary <- bind_rows(summary_rows)

pair_rows <- list(); j <- 0L
pairs <- combn(methods, 2, simplify = FALSE)
for (condition in c(conditions, "ALL")) for (pair in pairs) {
  d <- if (condition == "ALL") per_rep else per_rep[per_rep$condition == condition, ]
  a <- d[d$method == pair[[1]], c("rep", "condition", "abs_mean_error")]
  b <- d[d$method == pair[[2]], c("rep", "condition", "abs_mean_error")]
  z <- merge(a, b, by = c("rep", "condition"), suffixes = c("_a", "_b"))
  delta <- z$abs_mean_error_a - z$abs_mean_error_b
  wt <- suppressWarnings(wilcox.test(z$abs_mean_error_a, z$abs_mean_error_b,
                                     paired = TRUE, exact = FALSE))
  j <- j + 1L
  pair_rows[[j]] <- data.frame(condition = condition, method_a = pair[[1]], method_b = pair[[2]],
                               n_pairs = nrow(z), mean_delta_abs_mean_error = mean(delta),
                               median_delta_abs_mean_error = median(delta), p_value = wt$p.value,
                               stringsAsFactors = FALSE)
}
pairwise <- bind_rows(pair_rows)
pairwise$p_adjust_bh <- ave(pairwise$p_value, pairwise$condition,
                            FUN = function(x) p.adjust(x, method = "BH"))

write.table(manifest, file.path(out_dir, "prostate_three_arm_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(matched, file.path(out_dir, "prostate_three_arm_matched_events.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
saveRDS(matched, file.path(out_dir, "prostate_three_arm_matched_events.rds"), compress = "xz")
write.table(per_rep, file.path(out_dir, "prostate_three_arm_per_replicate.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_dir, "prostate_three_arm_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(pairwise, file.path(out_dir, "prostate_three_arm_pairwise_tests.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

plot_dat <- summary[summary$metric == "abs_mean_error", ]
pdf(file.path(out_dir, "prostate_three_arm_abs_mean_error.pdf"), width = 10, height = 5.5)
cols <- c("#31688e", "#35b779", "#d95f02"); names(cols) <- methods
matplot(seq_along(conditions), matrix(NA_real_, nrow = length(conditions), ncol = 3),
        type = "n", xaxt = "n", xlab = "Mixture condition", ylab = "Absolute mean CCF error",
        ylim = range(c(plot_dat$ci_lo, plot_dat$ci_hi), finite = TRUE))
axis(1, seq_along(conditions), conditions)
for (i in seq_along(methods)) {
  d <- plot_dat[plot_dat$method == methods[[i]], ]; d <- d[match(conditions, d$condition), ]
  xx <- seq_along(conditions) + (i - 2) * .08
  arrows(xx, d$ci_lo, xx, d$ci_hi, angle = 90, code = 3, length = .035, col = cols[[i]])
  lines(xx, d$estimate, type = "b", pch = 15 + i, col = cols[[i]])
}
legend("topright", methods, col = cols, pch = 16:18, lty = 1, bty = "n")
dev.off()

overall <- pairwise[pairwise$condition == "ALL", ]
overall_method <- aggregate(abs_mean_error ~ method, per_rep, mean)
report <- c(
  "# Prostate-mixture three-arm CCF comparison", "",
  sprintf("Generated: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
  "## Arms", "",
  "1. SVCFit.",
  "2. SVclone truth-CCF-assisted: known CCF selects multiplicity; this is not truth-CN assistance.",
  "3. SVclone non-assisted: native CCube endpoint CCF from the FACETS-based rerun.", "",
  "## Validation", "",
  sprintf("- 330/330 replicate-condition combinations passed."),
  sprintf("- %s event-condition observations were scored on the identical three-way shared set.",
          format(nrow(matched), big.mark = ",")),
  sprintf("- Shared events per replicate-condition: median %d, range %d-%d.",
          median(manifest$three_way_events), min(manifest$three_way_events),
          max(manifest$three_way_events)),
  "- Assisted and non-assisted SVclone purities agreed for every condition.", "",
  "## Overall absolute mean CCF error", "",
  "Each replicate-condition contributes equally.", "",
  "| Method | Mean absolute mean error |",
  "|---|---:|",
  vapply(seq_len(nrow(overall_method)), function(i) sprintf("| %s | %.5f |",
    overall_method$method[i], overall_method$abs_mean_error[i]), character(1)), "",
  "## Overall paired comparisons", "",
  "Differences are method A minus method B in replicate-condition absolute mean CCF error; negative favors method A.", "",
  paste0("| Method A | Method B | Pairs | Mean difference | BH-adjusted p |"),
  "|---|---|---:|---:|---:|",
  vapply(seq_len(nrow(overall)), function(i) sprintf("| %s | %s | %d | %.5f | %.4g |",
    overall$method_a[i], overall$method_b[i], overall$n_pairs[i],
    overall$mean_delta_abs_mean_error[i], overall$p_adjust_bh[i]), character(1)), "",
  "Condition-level estimates and bootstrap confidence intervals are in `prostate_three_arm_summary.tsv`.",
  "The analysis uses exact event identifiers, matching the established prostate benchmark's shared-event estimand.",
  "The three-way shared set requires every method to evaluate an event, so it is narrower than a two-way set; summary values from the two are not expected to be identical."
)
writeLines(report, file.path(out_dir, "PROSTATE_THREE_ARM_REPORT.md"))

git_commit <- tryCatch(system2("git", c("-C", shQuote(dirname(dirname(script_dir))), "rev-parse", "HEAD"),
                               stdout = TRUE, stderr = FALSE), error = function(e) "unknown")
prov <- c(sprintf("workflow_commit\t%s", git_commit[[1]]),
          sprintf("seed\t%d", seed), sprintf("bootstrap\t%d", bootstrap_n),
          sprintf("svcfit_chrx\t%s", normalizePath(opt[["svcfit-chrx"]])),
          sprintf("svcfit_45\t%s", normalizePath(opt[["svcfit-45"]])),
          sprintf("assisted_root\t%s", normalizePath(opt[["assisted-root"]])),
          sprintf("fair_root\t%s", normalizePath(opt[["fair-root"]])),
          sprintf("truth_dir\t%s", normalizePath(opt[["truth-dir"]])),
          capture.output(sessionInfo()))
writeLines(prov, file.path(out_dir, "PROVENANCE.txt"))
cat("Prostate three-arm analysis gate: PASS\nOutput: ", out_dir, "\n", sep = "")
