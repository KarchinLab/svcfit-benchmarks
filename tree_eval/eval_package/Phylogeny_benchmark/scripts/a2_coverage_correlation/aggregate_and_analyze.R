#!/usr/bin/env Rscript
#
# A2 coverage correlation — aggregation, statistics, and Fig A2-3.
#
# Reads all per-run CSVs from PER_RUN_DIR, joins with the existing
# signed-error table, then runs:
#   - Spearman rho stratified by clone, with bootstrap 95% CIs over
#     replicates (replicates resampled as units to preserve pre/post pairing).
#   - lmer(signed_error ~ log10(eff_cov) * clone + (1 | replicate)) with
#     interaction LRT against a no-interaction null.
# Writes Fig A2-3 (per-clone scatter), summary tables, and a joined dataset.

script_dir <- Sys.getenv("A2_SCRIPT_DIR", unset = ".")
source(file.path(script_dir, "config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Aggregate per-run outputs.
# -----------------------------------------------------------------------------

per_run_files <- list.files(PER_RUN_DIR, pattern = "^eff_cov_.*\\.csv$",
                            full.names = TRUE)
if (length(per_run_files) == 0L) stop("No per-run files found in ", PER_RUN_DIR)

eff_cov <- do.call(rbind,
  lapply(per_run_files, read.csv, stringsAsFactors = FALSE))
message(sprintf("Loaded %d (run x cluster) rows from %d per-run files.",
                nrow(eff_cov), length(per_run_files)))

# -----------------------------------------------------------------------------
# 2. Join with signed-error table.
# -----------------------------------------------------------------------------

if (!file.exists(SIGNED_ERROR_TABLE)) {
  stop("Signed-error table not found: ", SIGNED_ERROR_TABLE)
}
sgn <- read.csv(SIGNED_ERROR_TABLE, stringsAsFactors = FALSE)
required <- c("purity", "replicate", "timepoint", "true_clone",
              "inferred_cluster_id", "true_ccf", "inferred_ccf", "signed_error")
missing  <- setdiff(required, colnames(sgn))
if (length(missing) > 0L) {
  stop("signed-error table missing columns: ", paste(missing, collapse = ", "))
}

df <- inner_join(sgn, eff_cov,
                 by = c("purity", "replicate", "timepoint", "inferred_cluster_id"))
if (nrow(df) == 0L) stop("Join produced 0 rows. Check inferred_cluster_id matches between files.")

df$log10_eff_cov <- log10(df$mean_AO)
clone_levels <- c("Trunk", "Subclone 1", "Subclone 2")
df$clone <- factor(df$true_clone, levels = clone_levels)
if (any(is.na(df$clone))) {
  warning(sprintf("%d rows had a true_clone value not in {%s} — dropped.",
                  sum(is.na(df$clone)), paste(clone_levels, collapse = ", ")))
  df <- df[!is.na(df$clone), ]
}

joined_path <- file.path(OUT_DIR, "a2_eff_cov_signed_error_joined.csv")
write.csv(df, joined_path, row.names = FALSE)
message(sprintf("Joined dataset: %d rows -> %s", nrow(df), joined_path))

# -----------------------------------------------------------------------------
# 3. Stratified Spearman with bootstrap CIs over replicates.
# -----------------------------------------------------------------------------

set.seed(BOOTSTRAP_SEED)

bootstrap_spearman <- function(sub, B) {
  reps <- unique(sub$replicate)
  rho_obs <- suppressWarnings(
    cor(sub$mean_AO, sub$signed_error, method = "spearman")
  )
  boots <- vapply(seq_len(B), function(b) {
    sampled <- sample(reps, size = length(reps), replace = TRUE)
    rows <- do.call(rbind, lapply(sampled, function(r) sub[sub$replicate == r, ]))
    suppressWarnings(cor(rows$mean_AO, rows$signed_error, method = "spearman"))
  }, FUN.VALUE = numeric(1))
  ci <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
  list(rho = rho_obs, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]))
}

spearman_tbl <- do.call(rbind, lapply(clone_levels, function(cl) {
  sub <- df[df$clone == cl, ]
  if (nrow(sub) < 3L) {
    return(data.frame(clone = cl, rho = NA, ci_lo = NA, ci_hi = NA, n = nrow(sub)))
  }
  res <- bootstrap_spearman(sub, BOOTSTRAP_B)
  data.frame(clone = cl, rho = res$rho, ci_lo = res$ci_lo, ci_hi = res$ci_hi,
             n = nrow(sub))
}))
spearman_path <- file.path(OUT_DIR, "a2_spearman_by_clone.csv")
write.csv(spearman_tbl, spearman_path, row.names = FALSE)
message("Stratified Spearman:")
print(spearman_tbl)

# -----------------------------------------------------------------------------
# 4. Mixed-effects regression with interaction LRT.
# -----------------------------------------------------------------------------

m_full <- lmer(signed_error ~ log10_eff_cov * clone + (1 | replicate), data = df)
m_null <- lmer(signed_error ~ log10_eff_cov + clone + (1 | replicate), data = df)
lrt    <- anova(m_null, m_full)

lmer_path <- file.path(OUT_DIR, "a2_lmer_summary.txt")
sink(lmer_path)
cat("Full model: signed_error ~ log10(eff_cov) * clone + (1 | replicate)\n\n")
print(summary(m_full))
cat("\n--------------------------------------------------------------\n")
cat("Null model: signed_error ~ log10(eff_cov) + clone + (1 | replicate)\n\n")
print(summary(m_null))
cat("\n--------------------------------------------------------------\n")
cat("Likelihood-ratio test (interaction):\n\n")
print(lrt)
sink()

interaction_p <- lrt$`Pr(>Chisq)`[2]
message(sprintf("Mixed model: interaction LRT p = %s -> %s",
                signif(interaction_p, 3), lmer_path))

# -----------------------------------------------------------------------------
# 5. Figure A2-3.
# -----------------------------------------------------------------------------

clone_colors <- c("Trunk" = "grey20",
                  "Subclone 1" = "steelblue",
                  "Subclone 2" = "firebrick")

p <- ggplot(df, aes(x = log10_eff_cov, y = signed_error, color = clone)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  scale_color_manual(values = clone_colors) +
  labs(x = expression(log[10] ~ "effective per-clone supporting coverage (mean AO)"),
       y = "Signed CCF error (true - inferred)",
       color = "Clone",
       title = "Figure A2-3. Per-clone signed CCF error vs. coverage") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(size = 12, face = "bold"))

fig_path <- file.path(OUT_DIR, "fig_a2_3.pdf")
ggsave(fig_path, p, width = 7, height = 5)
message("Wrote ", fig_path)

# -----------------------------------------------------------------------------
# 6. One-line interpretation against the protocol's decision rules.
# -----------------------------------------------------------------------------

interpret <- function(spearman_tbl, p_int) {
  rho_t  <- spearman_tbl$rho[spearman_tbl$clone == "Trunk"]
  ci_t   <- c(spearman_tbl$ci_lo[spearman_tbl$clone == "Trunk"],
              spearman_tbl$ci_hi[spearman_tbl$clone == "Trunk"])
  rho_s1 <- spearman_tbl$rho[spearman_tbl$clone == "Subclone 1"]
  ci_s1  <- c(spearman_tbl$ci_lo[spearman_tbl$clone == "Subclone 1"],
              spearman_tbl$ci_hi[spearman_tbl$clone == "Subclone 1"])
  rho_s2 <- spearman_tbl$rho[spearman_tbl$clone == "Subclone 2"]
  ci_s2  <- c(spearman_tbl$ci_lo[spearman_tbl$clone == "Subclone 2"],
              spearman_tbl$ci_hi[spearman_tbl$clone == "Subclone 2"])

  trunk_flat   <- !is.na(rho_t) && rho_t > -0.2 && rho_t < 0.2 &&
                  ci_t[1] < 0 && ci_t[2] > 0
  s1_negative  <- !is.na(rho_s1) && rho_s1 < -0.3 && ci_s1[2] < 0
  s2_negative  <- !is.na(rho_s2) && rho_s2 < -0.3 && ci_s2[2] < 0
  interaction_sig <- !is.na(p_int) && p_int < 0.05

  if (trunk_flat && s1_negative && s2_negative && interaction_sig) {
    "CONFIRMED: dual-mechanism (subclones coverage-limited; trunk one-sided censored)."
  } else if (s1_negative && s2_negative && !trunk_flat) {
    "PARTIAL: subclones coverage-limited as predicted, but trunk also shows non-flat correlation -- reconsider trunk mechanism (SV-type-specific bias?)."
  } else if (!s1_negative || !s2_negative) {
    "WEAK: subclone coverage relationship weaker than expected -- consider switching to definition (B) breakpoint depth, or check SV-type composition."
  } else {
    "AMBIGUOUS: pattern does not match any pre-specified decision rule cleanly. Inspect Fig A2-3 manually."
  }
}

verdict <- interpret(spearman_tbl, interaction_p)
verdict_path <- file.path(OUT_DIR, "a2_interpretation.txt")
writeLines(verdict, verdict_path)
message(sprintf("Interpretation: %s\n  -> %s", verdict, verdict_path))
