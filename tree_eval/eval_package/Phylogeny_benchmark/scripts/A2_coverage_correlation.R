#!/usr/bin/env Rscript
###############################################################################
# A2_coverage_correlation.R
#
# Tests two mechanisms behind the purity-dependent CCF error pattern (Fig A2-1B):
#   (1) Coverage-limited inference for subclones (predicts strong negative rho)
#   (2) One-sided censoring at CCF=1 for the trunk (predicts near-zero rho)
#
# Effective per-clone coverage = mean SVtyper AO across SVs assigned to that
# clone (definition A from the protocol; uses inferred cluster assignment).
#
# Outputs:
#   coverage_correlation_data.csv   Long-format data (one row per rep×clone×tp)
#   spearman_results.csv            rho + 95% bootstrap CI per clone
#   mixed_model_coefs.csv           lmer coefficient table
#   fig_a2_3.pdf                    Figure A2-3
#
# Usage:
#   Rscript A2_coverage_correlation.R \
#     --eval_csv  /path/to/evaluation_summary.csv \
#     --work_dir  "$TREE_EVAL_LONGITUDINAL" \
#     --out_dir   /path/to/output
###############################################################################

library(optparse)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lme4)

option_list <- list(
  make_option(c("--eval_csv"),   type = "character",
              help = "Path to evaluation_summary.csv"),
  make_option(c("--work_dir"),   type = "character",
              help = "Longitudinal simulation base dir (parent of S1/)"),
  make_option(c("--out_dir"),    type = "character",
              help = "Output directory for all results"),
  make_option(c("--cov"),        type = "integer",   default = 50L,
              help = "Nominal coverage used in simulation [default: 50]"),
  make_option(c("--purities"),   type = "character", default = "20,40,60,80",
              help = "Comma-separated purities to include [default: 20,40,60,80]"),
  make_option(c("--n_iter"),     type = "integer",   default = 1000L,
              help = "Bootstrap iterations for 95% CI [default: 1000]"),
  make_option(c("--pre_stage"),  type = "character", default = "pre_BAT",
              help = "Stage label for pre-treatment in sv2cluster.csv"),
  make_option(c("--post_stage"), type = "character", default = "on_BAT",
              help = "Stage label for post-treatment in sv2cluster.csv"),
  make_option(c("--seed"),       type = "integer",   default = 42L,
              help = "Random seed for bootstrap [default: 42]")
)

parser <- OptionParser(option_list = option_list, add_help_option = TRUE)
opts   <- parse_args(parser)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(opts$seed)

purities <- as.integer(strsplit(opts$purities, ",")[[1]])

###############################################################################
# Ground-truth CCFs: replicate bash integer division from longi_short.sh
###############################################################################
true_clone_fracs_s1 <- function(pur) {
  idiv <- function(a, b) as.integer(a) %/% as.integer(b)
  s1_t1_c2 <- idiv(pur * 5L, 6L); s1_t1_c3 <- pur - s1_t1_c2
  s1_t2_c2 <- idiv(pur * 1L, 6L); s1_t2_c3 <- pur - s1_t2_c2
  list(
    t1 = c(trunk = 1.0, sub1 = s1_t1_c2 / pur, sub2 = s1_t1_c3 / pur),
    t2 = c(trunk = 1.0, sub1 = s1_t2_c2 / pur, sub2 = s1_t2_c3 / pur)
  )
}

###############################################################################
# Parse SVtyper VCF and extract AO from the first bulk column.
#
# The per-SV genotyped VCFs carry two identical 'bulk' columns; the first is
# used. AO is the SVtyper alternate-allele observation count (split + PE +
# clipped-read evidence, post-genotyping). Manta's PR/SR are NOT used.
#
# Returns: data.frame(sv_id, AO)
###############################################################################
parse_vcf_ao <- function(vcf_file) {
  lines       <- readLines(vcf_file, warn = FALSE)
  header_line <- lines[startsWith(lines, "#CHROM")][1]
  data_lines  <- lines[!startsWith(lines, "#")]

  if (length(data_lines) == 0 || is.na(header_line))
    return(data.frame(sv_id = character(0), AO = integer(0),
                      stringsAsFactors = FALSE))

  col_names  <- strsplit(header_line, "\t")[[1]]
  fmt_idx    <- which(col_names == "FORMAT")
  samp_idx   <- fmt_idx + 1L    # first 'bulk' column

  rows  <- strsplit(data_lines, "\t")
  sv_id <- vapply(rows, `[[`, character(1), 3L)       # VCF ID column
  fmt   <- strsplit(vapply(rows, `[[`, character(1), fmt_idx),  ":")
  samp  <- strsplit(vapply(rows, `[[`, character(1), samp_idx), ":")

  ao_val <- mapply(function(f, s) {
    idx <- match("AO", f)
    if (is.na(idx) || idx > length(s)) return(NA_integer_)
    suppressWarnings(as.integer(s[[idx]]))
  }, fmt, samp, SIMPLIFY = TRUE)

  data.frame(sv_id = sv_id, AO = ao_val, stringsAsFactors = FALSE)
}

###############################################################################
# Assign sub1/sub2 labels to non-trunk clusters using the same Hungarian
# matching as evaluate_downstream.R (minimises total CCF MAE).
#
# centroids: data.frame(cluster_num, t1_ccf, t2_ccf)
# trunk_cluster: string like "cluster3"
# true_ccfs: list(t1=named numeric, t2=named numeric) from true_clone_fracs_s1
#
# Returns: data.frame(cluster_num, clone)   [3 rows when successful]
###############################################################################
assign_clone_labels <- function(centroids, trunk_cluster, true_ccfs) {
  non_trunk <- centroids[centroids$cluster_num != trunk_cluster, , drop = FALSE]

  trunk_row <- data.frame(cluster_num = trunk_cluster, clone = "trunk",
                          stringsAsFactors = FALSE)
  if (nrow(non_trunk) != 2L)
    return(trunk_row)    # fallback: only trunk labelled; sub1/sub2 unresolved

  cands <- list(
    list(sub1 = non_trunk$cluster_num[1L], sub2 = non_trunk$cluster_num[2L]),
    list(sub1 = non_trunk$cluster_num[2L], sub2 = non_trunk$cluster_num[1L])
  )

  best_err  <- Inf
  best_cand <- cands[[1]]

  for (cand in cands) {
    inf_t1 <- c(
      trunk = centroids$t1_ccf[centroids$cluster_num == trunk_cluster],
      sub1  = centroids$t1_ccf[centroids$cluster_num == cand$sub1],
      sub2  = centroids$t1_ccf[centroids$cluster_num == cand$sub2]
    )
    inf_t2 <- c(
      trunk = centroids$t2_ccf[centroids$cluster_num == trunk_cluster],
      sub1  = centroids$t2_ccf[centroids$cluster_num == cand$sub1],
      sub2  = centroids$t2_ccf[centroids$cluster_num == cand$sub2]
    )
    err <- sum(abs(inf_t1 - true_ccfs$t1)) + sum(abs(inf_t2 - true_ccfs$t2))
    if (err < best_err) {
      best_err  <- err
      best_cand <- cand
    }
  }

  bind_rows(trunk_row,
    data.frame(
      cluster_num = c(best_cand$sub1,  best_cand$sub2),
      clone       = c("sub1",          "sub2"),
      stringsAsFactors = FALSE
    )
  )
}

###############################################################################
# Compute effective coverage for one replicate (boot × purity × experiment).
#
# Returns a data.frame with columns:
#   clone, timepoint, eff_cov, n_svs
# or NULL on any unrecoverable failure (missing files, wrong cluster count).
###############################################################################
compute_eff_cov <- function(boot, purity, exp_name,
                             work_dir, cov,
                             pre_stage_lbl, post_stage_lbl,
                             true_ccfs) {
  longi_dir <- if (boot == 0L) {
    file.path(work_dir, "S1", paste0("c", cov, "p", purity), exp_name)
  } else {
    file.path(work_dir, "S1", paste0("c", cov, "p", purity),
              paste0("b", boot), exp_name)
  }

  cluster_file   <- file.path(longi_dir, "clustering", "sv2cluster.csv")
  centroids_file <- file.path(longi_dir, "clustering", "cluster_centroids.csv")
  tree_file      <- file.path(longi_dir, "tree",       "tree_result.rds")

  if (!file.exists(cluster_file) || !file.exists(centroids_file) ||
      !file.exists(tree_file))
    return(NULL)

  sv2cluster <- read_csv(cluster_file,   show_col_types = FALSE)
  centroids  <- read_csv(centroids_file, show_col_types = FALSE) %>%
    filter(pair == 1) %>%
    select(cluster_num, t1_ccf = f_pre, t2_ccf = f_day85)
  tree_res   <- readRDS(tree_file)
  tree_edges <- tree_res[[1]]

  trunk_idx <- tree_edges$child[tree_edges$parent == "root"]
  if (length(trunk_idx) != 1L) return(NULL)
  trunk_cluster <- paste0("cluster", trunk_idx)

  clone_map <- assign_clone_labels(centroids, trunk_cluster, true_ccfs)
  if (nrow(clone_map) < 3L) return(NULL)   # sub1/sub2 unresolved (e.g. 2-cluster solution)

  tp_specs <- list(
    list(stage_lbl = pre_stage_lbl,  tp_str = "t1", timepoint = "pre"),
    list(stage_lbl = post_stage_lbl, tp_str = "t2", timepoint = "post")
  )

  result_rows <- lapply(tp_specs, function(spec) {
    vcf_file <- file.path(longi_dir, spec$tp_str, "svtyp",
                          sprintf("svt_%s_%s.vcf", exp_name, spec$tp_str))
    if (!file.exists(vcf_file)) return(NULL)

    vcf_ao <- parse_vcf_ao(vcf_file)
    if (nrow(vcf_ao) == 0L) return(NULL)

    sv_tp <- sv2cluster %>%
      filter(stage == spec$stage_lbl) %>%
      left_join(clone_map, by = "cluster_num") %>%
      filter(!is.na(clone)) %>%
      left_join(vcf_ao, by = c("shared" = "sv_id"))

    sv_tp %>%
      group_by(clone) %>%
      summarise(
        eff_cov  = mean(AO, na.rm = TRUE),
        n_svs    = sum(!is.na(AO)),
        .groups  = "drop"
      ) %>%
      mutate(timepoint = spec$timepoint)
  })

  bind_rows(result_rows[!vapply(result_rows, is.null, logical(1))])
}

###############################################################################
# Main: iterate over replicates, collect coverage + signed errors
###############################################################################
cat("Reading evaluation summary...\n")
eval_sum <- read_csv(opts$eval_csv, show_col_types = FALSE) %>%
  filter(purity %in% purities)

cat(sprintf("Replicates to process: %d\n", nrow(eval_sum)))

cov_rows <- list()

for (i in seq_len(nrow(eval_sum))) {
  row      <- eval_sum[i, ]
  boot_i   <- row$boot
  pur_i    <- row$purity
  exp_i    <- row$experiment
  true_ccfs <- true_clone_fracs_s1(pur_i)

  cov_df <- compute_eff_cov(
    boot_i, pur_i, exp_i,
    opts$work_dir, opts$cov,
    opts$pre_stage, opts$post_stage,
    true_ccfs
  )

  if (is.null(cov_df) || nrow(cov_df) == 0L) {
    cat(sprintf("  Skipped boot=%d p=%d %s\n", boot_i, pur_i, exp_i))
    next
  }

  cov_df <- cov_df %>%
    mutate(
      boot       = boot_i,
      purity     = pur_i,
      experiment = exp_i,
      replicate  = paste0("b", boot_i, "_p", pur_i, "_", exp_i)
    )
  cov_rows <- c(cov_rows, list(cov_df))
}

if (length(cov_rows) == 0L) {
  cat("No coverage data computed. Exiting.\n")
  quit(status = 1)
}

cov_all <- bind_rows(cov_rows)

###############################################################################
# Reshape evaluation_summary signed errors to long format and join
###############################################################################
err_long <- eval_sum %>%
  select(boot, purity, experiment, starts_with("ccf_err_")) %>%
  pivot_longer(
    cols       = starts_with("ccf_err_"),
    names_to   = "metric",
    values_to  = "signed_error"
  ) %>%
  mutate(
    clone     = sub("ccf_err_(trunk|sub1|sub2)_(pre|post)", "\\1", metric),
    timepoint = sub("ccf_err_(trunk|sub1|sub2)_(pre|post)", "\\2", metric)
  ) %>%
  select(boot, purity, experiment, clone, timepoint, signed_error)

df <- cov_all %>%
  left_join(err_long,
            by = c("boot", "purity", "experiment", "clone", "timepoint")) %>%
  filter(!is.na(signed_error), !is.na(eff_cov), eff_cov > 0, n_svs > 0) %>%
  mutate(
    log10_eff_cov = log10(eff_cov),
    clone = factor(clone,
                   levels = c("trunk", "sub1", "sub2"),
                   labels = c("Trunk", "Sub1", "Sub2")),
    timepoint = factor(timepoint,
                       levels = c("pre", "post"),
                       labels = c("Pre", "Post"))
  )

cat(sprintf("\nBuilt long-format data frame: %d rows (expected up to 480)\n", nrow(df)))
write_csv(df, file.path(opts$out_dir, "coverage_correlation_data.csv"))

###############################################################################
# Spearman rho per clone — bootstrap 95% CI resampling over replicate IDs.
# A "replicate" = one (boot, purity, experiment) triple contributing exactly
# 2 obs per clone (Pre + Post). The paired Pre/Post structure is preserved
# by resampling whole replicate blocks.
###############################################################################
cat("\n--- Spearman rho (stratified by clone) ---\n")

n_iter <- opts$n_iter

spearman_results <- lapply(levels(df$clone), function(cl) {
  sub_df      <- df %>% filter(clone == cl)
  rep_ids     <- unique(sub_df$replicate)

  obs_rho <- cor(sub_df$log10_eff_cov, sub_df$signed_error,
                 method = "spearman", use = "complete.obs")

  boot_rhos <- replicate(n_iter, {
    sampled <- sample(rep_ids, length(rep_ids), replace = TRUE)
    bs_df   <- bind_rows(lapply(sampled,
                                function(r) sub_df[sub_df$replicate == r, ]))
    cor(bs_df$log10_eff_cov, bs_df$signed_error,
        method = "spearman", use = "complete.obs")
  })

  ci <- quantile(boot_rhos, c(0.025, 0.975), na.rm = TRUE)

  data.frame(
    clone   = cl,
    rho     = obs_rho,
    ci_lo   = ci[[1]],
    ci_hi   = ci[[2]],
    n_obs   = nrow(sub_df),
    n_reps  = length(rep_ids),
    stringsAsFactors = FALSE
  )
})

spearman_df <- bind_rows(spearman_results)
cat("\nSpearman rho [95% bootstrap CI]:\n")
print(spearman_df, digits = 3)
write_csv(spearman_df, file.path(opts$out_dir, "spearman_results.csv"))

###############################################################################
# Mixed-effects regression (required)
#
# Model: signed_error ~ log10(eff_cov) * clone + (1 | replicate)
# Null:  signed_error ~ log10(eff_cov) + clone  + (1 | replicate)
# LRT of interaction tests whether the coverage-error relationship is
# clone-specific (the headline claim).
###############################################################################
cat("\n--- Mixed-effects regression ---\n")
cat("Full:  signed_error ~ log10_eff_cov * clone + (1 | replicate)\n")
cat("Null:  signed_error ~ log10_eff_cov + clone  + (1 | replicate)\n\n")

m_full <- lmer(
  signed_error ~ log10_eff_cov * clone + (1 | replicate),
  data = df, REML = FALSE
)
m_null <- lmer(
  signed_error ~ log10_eff_cov + clone + (1 | replicate),
  data = df, REML = FALSE
)

lrt      <- anova(m_null, m_full)
lrt_pval <- lrt[["Pr(>Chisq)"]][2]

cat("Full model summary:\n")
print(summary(m_full))
cat("\nLikelihood-ratio test (interaction vs. no interaction):\n")
print(lrt)
cat(sprintf("\nInteraction LRT p-value: %.4g\n", lrt_pval))

# Derive per-clone slopes from fixed effects
coefs     <- as.data.frame(coef(summary(m_full)))
coefs$term <- rownames(coefs)
write_csv(coefs, file.path(opts$out_dir, "mixed_model_coefs.csv"))

get_slope <- function(coefs, clone_label) {
  base_slope  <- coefs$Estimate[coefs$term == "log10_eff_cov"]
  inter_term  <- paste0("log10_eff_cov:clone", clone_label)
  interaction <- if (inter_term %in% coefs$term)
    coefs$Estimate[coefs$term == inter_term] else 0
  base_slope + interaction
}

cat("\nPer-clone slopes of log10(eff_cov):\n")
for (cl in levels(df$clone)) {
  cat(sprintf("  %s: %.4f\n", cl, get_slope(coefs, cl)))
}

###############################################################################
# Figure A2-3: signed CCF error vs. log10(eff_cov), per-clone LOESS
###############################################################################
clone_colors <- c("Trunk" = "#4d4d4d", "Sub1" = "#1f78b4", "Sub2" = "#e31a1c")

rho_labels <- spearman_df %>%
  mutate(
    clone = factor(clone, levels = c("Trunk", "Sub1", "Sub2")),
    label = sprintf("rho == %.2f~'['*%.2f*','~%.2f*']'", rho, ci_lo, ci_hi)
  )

# x position for annotation: 5% from left of each panel's range
label_x <- df %>%
  group_by(clone) %>%
  summarise(x_pos = min(log10_eff_cov, na.rm = TRUE), .groups = "drop")

rho_labels <- rho_labels %>%
  left_join(label_x, by = "clone")

p <- df %>%
  ggplot(aes(x = log10_eff_cov, y = signed_error,
             color = clone, fill = clone)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(alpha = 0.40, size = 1.6, shape = 16) +
  geom_smooth(method = "loess", formula = y ~ x,
              se = TRUE, alpha = 0.15, linewidth = 0.9) +
  geom_text(
    data    = rho_labels,
    aes(x = x_pos, y = Inf, label = label),
    hjust = 0, vjust = 1.4, size = 3.2, color = "black",
    parse = TRUE, inherit.aes = FALSE
  ) +
  scale_color_manual(values = clone_colors) +
  scale_fill_manual( values = clone_colors) +
  facet_wrap(~clone, scales = "free_x") +
  labs(
    x       = expression(log[10](mean~AO~per~clone)),
    y       = "Signed CCF error (true − inferred)",
    title   = "Figure A2-3. Per-clone signed CCF error vs. effective coverage",
    caption = sprintf(
      paste0(
        "Each point is one (replicate × clone × timepoint) observation, pooled across purities and Pre/Post timepoints.\n",
        "Lines are per-clone LOESS fits (shaded = 95%% CI). ρ = Spearman rho [95%% bootstrap CI, %d iterations].\n",
        "n = %d observations total."
      ),
      n_iter, nrow(df)
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 12),
    plot.caption     = element_text(hjust = 0, size = 8)
  )

ggsave(file.path(opts$out_dir, "fig_a2_3.pdf"), p, width = 10, height = 4.5)
cat("\nFigure A2-3 saved.\n")

###############################################################################
# Decision rules
###############################################################################
cat("\n--- Decision rules ---\n")

for (cl in c("Sub1", "Sub2", "Trunk")) {
  r <- spearman_df[spearman_df$clone == cl, ]
  if (nrow(r) == 0L) { cat(sprintf("  %s: no data\n", cl)); next }

  if (cl %in% c("Sub1", "Sub2")) {
    met <- r$rho < -0.3 && r$ci_hi < 0
    cat(sprintf(
      "  %s: rho = %+.3f [%+.3f, %+.3f]  | predicted rho < -0.3 and CI excludes 0: %s\n",
      cl, r$rho, r$ci_lo, r$ci_hi, if (met) "YES" else "NO"
    ))
  } else {
    met <- abs(r$rho) < 0.2 && r$ci_lo < 0 && r$ci_hi > 0
    cat(sprintf(
      "  Trunk: rho = %+.3f [%+.3f, %+.3f]  | predicted |rho| < 0.2 and CI includes 0: %s\n",
      r$rho, r$ci_lo, r$ci_hi, if (met) "YES" else "NO"
    ))
  }
}

cat(sprintf(
  "  Interaction LRT: p = %.4g  | predicted p < 0.05: %s\n",
  lrt_pval, if (lrt_pval < 0.05) "YES" else "NO"
))

overall <- {
  sub1  <- spearman_df[spearman_df$clone == "Sub1",  ]
  sub2  <- spearman_df[spearman_df$clone == "Sub2",  ]
  trunk <- spearman_df[spearman_df$clone == "Trunk", ]
  nrow(sub1) > 0 && nrow(sub2) > 0 && nrow(trunk) > 0 &&
    sub1$rho < -0.3 && sub1$ci_hi < 0 &&
    sub2$rho < -0.3 && sub2$ci_hi < 0 &&
    abs(trunk$rho) < 0.2 && trunk$ci_lo < 0 && trunk$ci_hi > 0 &&
    lrt_pval < 0.05
}
cat(sprintf("\nAll four decision rules met (dual-mechanism confirmed): %s\n",
            if (overall) "YES" else "NO"))

cat(sprintf("\nAll outputs written to: %s\n", opts$out_dir))
cat("Done.\n")
