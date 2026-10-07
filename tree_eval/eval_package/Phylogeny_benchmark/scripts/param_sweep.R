#!/usr/bin/env Rscript
###############################################################################
# param_sweep.R
#
# Two modes:
#   --task <array_id>   Run one SLURM array task (called by param_sweep_slurm.sh)
#   --merge             Collect summary_row.csv files into summary.csv
#
# Submit: sbatch param_sweep_slurm.sh
# Merge:  Rscript param_sweep.R --merge
###############################################################################

# ── CONFIG ────────────────────────────────────────────────────────────────────
.args        <- commandArgs(trailingOnly = FALSE)
.script_file <- sub("--file=", "", .args[grep("--file=", .args)])
.script_dir  <- if (length(.script_file) > 0) dirname(normalizePath(.script_file)) else getwd()
.pub_dir     <- dirname(.script_dir)

WORK_DIR   <- file.path(.pub_dir, "output")
SWEEP_DIR  <- file.path(WORK_DIR, "param_sweep")
PYTHON_BIN <- Sys.getenv("RETICULATE_PYTHON")
if (!nzchar(PYTHON_BIN) || !file.exists(PYTHON_BIN))
  stop("RETICULATE_PYTHON must name the configured project Python executable.")
COV        <- 50L

CONCENTRATION_VALS <- c(0.1, 0.5, 1, 3, 5, 10)
MIN_DIST_VALS      <- c(0.1, 0.2, 0.3)

TEST_CASES <- list(
  list(scenario = "S3", purity = 40L, exp = "exp2", note = "1 cluster — extreme under-clustering"),
  list(scenario = "S3", purity = 40L, exp = "exp1", note = "2 clusters — sub-clones merged"),
  list(scenario = "S3", purity = 60L, exp = "exp2", note = "3 clusters but linear chain topology"),
  list(scenario = "S3", purity = 80L, exp = "exp2", note = "3 clusters — wrong topology at high purity"),
  list(scenario = "S1", purity = 40L, exp = "exp2", note = "3 clusters but linear chain topology"),
  list(scenario = "S1", purity = 60L, exp = "exp2", note = "4 clusters — over-clustering"),
  list(scenario = "S1", purity = 80L, exp = "exp2", note = "4 clusters — over-clustering"),
  list(scenario = "S4", purity = 40L, exp = "exp2", note = "2 clusters — sub-clones merged"),
  list(scenario = "S4", purity = 60L, exp = "exp2", note = "2 clusters — sub-clones merged"),
  list(scenario = "S4", purity = 80L, exp = "exp2", note = "4 clusters — over-clustering")
)

# ── PACKAGES (always) ─────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ── HELPERS ───────────────────────────────────────────────────────────────────
is_correct_topology <- function(tree_edges) {
  if (is.null(tree_edges) || nrow(tree_edges) == 0) return(FALSE)
  root_children <- tree_edges$child[tree_edges$parent == "root"]
  if (length(root_children) != 1) return(FALSE)
  trunk <- root_children[1]
  trunk_children <- tree_edges$child[tree_edges$parent == trunk]
  if (length(trunk_children) != 2) return(FALSE)
  all(vapply(trunk_children,
             function(n) length(tree_edges$child[tree_edges$parent == n]) == 0,
             logical(1)))
}

topology_string <- function(tree_edges) {
  if (is.null(tree_edges) || nrow(tree_edges) == 0) return(NA_character_)
  paste(apply(tree_edges, 1, function(r) paste(r["parent"], r["child"], sep = "->")),
        collapse = "; ")
}

make_scatter <- function(cluster_result, scenario, purity, exp, conc, md) {
  centroids_df <- cluster_result %>% distinct(cluster_num, pre_center, post_center)
  ggplot(cluster_result, aes(x = pre_ccf, y = post_ccf, colour = cluster_num)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_point(data = centroids_df,
               aes(x = pre_center, y = post_center, colour = cluster_num),
               shape = 8, size = 5, stroke = 1.5, show.legend = FALSE) +
    geom_text(data = centroids_df,
              aes(x = pre_center, y = post_center,
                  label = sprintf("(%0.3f, %0.3f)", pre_center, post_center)),
              colour = "black", size = 3, vjust = -1.2, hjust = 0.5) +
    scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    labs(
      title   = sprintf("%s / p%s / %s", scenario, purity, exp),
      x       = "t1 CCF", y = "t2 CCF", colour = "Cluster",
      caption = sprintf("concentration = %s  |  min_dist = %s  |  n_clusters = %d",
                        conc, md, nrow(centroids_df))
    ) +
    theme_bw(base_size = 12) +
    theme(plot.caption = element_text(hjust = 0, size = 9))
}

make_delta_scatter <- function(cluster_result, scenario, purity, exp, conc, md) {
  df <- cluster_result %>% mutate(delta_ccf = post_ccf - pre_ccf)
  centroids_df <- cluster_result %>%
    distinct(cluster_num, pre_center, post_center) %>%
    mutate(delta_center = post_center - pre_center)
  ggplot(df, aes(x = delta_ccf, y = pre_ccf, colour = cluster_num)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_point(data = centroids_df,
               aes(x = delta_center, y = pre_center, colour = cluster_num),
               shape = 8, size = 5, stroke = 1.5, show.legend = FALSE) +
    geom_text(data = centroids_df,
              aes(x = delta_center, y = pre_center,
                  label = sprintf("(%0.3f, %0.3f)", delta_center, pre_center)),
              colour = "black", size = 3, vjust = -1.2, hjust = 0.5) +
    scale_x_continuous(limits = c(-1.05, 1.05), breaks = seq(-1, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    labs(
      title   = sprintf("%s / p%s / %s", scenario, purity, exp),
      x       = "t2 − t1 CCF", y = "t1 CCF", colour = "Cluster",
      caption = sprintf("concentration = %s  |  min_dist = %s  |  n_clusters = %d",
                        conc, md, nrow(centroids_df))
    ) +
    theme_bw(base_size = 12) +
    theme(plot.caption = element_text(hjust = 0, size = 9))
}

# ── TASK ──────────────────────────────────────────────────────────────────────
run_task <- function(task_id) {
  suppressPackageStartupMessages({
    library(tidyr)
    library(ggplot2)
    library(reticulate)
    library(SVCFit)
  })

  # Map task_id → (case, conc, md)
  n_combos <- length(CONCENTRATION_VALS) * length(MIN_DIST_VALS)  # 40
  case_idx <- task_id %/% n_combos + 1L
  combo_idx <- task_id %% n_combos
  conc_idx <- combo_idx %/% length(MIN_DIST_VALS) + 1L
  md_idx   <- combo_idx %%  length(MIN_DIST_VALS) + 1L

  tc       <- TEST_CASES[[case_idx]]
  scenario <- tc$scenario
  purity   <- tc$purity
  exp_name <- tc$exp
  note     <- tc$note
  conc     <- CONCENTRATION_VALS[conc_idx]
  md       <- MIN_DIST_VALS[md_idx]

  message(sprintf("Task %d: %s / p%s / %s  |  concentration=%.1f  min_dist=%.2f",
                  task_id, scenario, purity, exp_name, conc, md))

  use_python(PYTHON_BIN, required = TRUE)

  case_label  <- sprintf("%s_p%s_%s", scenario, purity, exp_name)
  combo_label <- sprintf("conc%.1f_md%.2f", conc, md)
  out_dir     <- file.path(SWEEP_DIR, case_label, combo_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  svcfit_dir <- file.path(SWEEP_DIR, case_label, "svcfit_output")
  combat_dir <- file.path(svcfit_dir, "COMBAT", "SVCFit_output")

  pur_frac  <- purity / 100
  purity_t2 <- if (scenario == "S4") pur_frac / 2 else pur_frac

  # ── Stage 1: SVCFit ──────────────────────────────────────────────────────────
  # combo_idx == 0 runs SVCFit; the other 17 combos per case wait for it.
  t1_bed <- file.path(svcfit_dir, sprintf("%s_p%s_%s_t1.bed", scenario, purity, exp_name))
  if (combo_idx == 0L) {
    if (!file.exists(t1_bed)) {
      message("  Running SVCFit ...")
      dir.create(combat_dir, recursive = TRUE, showWarnings = FALSE)
      longi_dir <- file.path(WORK_DIR, scenario, paste0("c", COV, "p", purity), exp_name)
      tp_input <- function(tp) list(
        het_snp = file.path(longi_dir, tp, "SNP",   sprintf("het_near_sv_%s_%s.vcf", exp_name, tp)),
        on_sv   = file.path(longi_dir, tp, "SNP",   sprintf("het_on_sv_%s_%s.vcf",  exp_name, tp)),
        sv      = file.path(longi_dir, tp, "svtyp", sprintf("svt_%s_%s.vcf",        exp_name, tp)),
        cnv     = file.path(longi_dir, tp, "facet", sprintf("%s_%s.bed",             exp_name, tp))
      )
      for (tp in c("t1", "t2")) {
        sn  <- sprintf("%s_p%s_%s_%s", scenario, purity, exp_name, tp)
        p   <- tp_input(tp)
        res <- run_svcfit(p_het = p$het_snp, p_onsv = p$on_sv, p_sv = p$sv, p_cnv = p$cnv,
                          samp = sn, exper = exp_name, thresh = 0.1)
        readr::write_delim(res$svcf, file.path(combat_dir, sprintf("%s.bed", sn)),
                           delim = "\t", quote = "none")
        readr::write_delim(res$svcf, file.path(svcfit_dir, sprintf("%s.bed", sn)),
                           delim = "\t", quote = "none")
        message(sprintf("    Written %d SVs: %s", nrow(res$svcf), sn))
      }
    } else {
      message("  SVCFit output already exists, skipping.")
    }
  } else {
    deadline <- Sys.time() + 30 * 60
    while (!file.exists(t1_bed)) {
      if (Sys.time() > deadline) stop("Timed out waiting for SVCFit output: ", t1_bed)
      Sys.sleep(30)
    }
  }

  pair_path <- file.path(out_dir, "pair_path.txt")
  write.table(
    data.frame(pre_BAT = sprintf("%s_p%s_%s_t1", scenario, purity, exp_name),
               on_BAT  = sprintf("%s_p%s_%s_t2", scenario, purity, exp_name)),
    pair_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  pur_path <- file.path(out_dir, "pur_path.txt")
  write.table(
    data.frame(sample = sprintf("%s_p%s_%s_%s", scenario, purity, exp_name, c("t1", "t2")),
               purity = c(pur_frac, purity_t2)),
    pur_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  result <- tryCatch({
    build_trees(
      run_clustering            = TRUE,
      pair_path                 = pair_path,
      pur_path                  = pur_path,
      data_dir                  = svcfit_dir,
      pair_num                  = 1L,
      deduplicate               = TRUE,
      concentration             = conc,
      min_dist                  = md,
      run_tree                  = TRUE,
      lineage_precedence_thresh = 0.2,
      sum_filter_thresh         = 0.35,
      ccf_floor                 = 0.1
    )
  }, error = function(e) list(error = conditionMessage(e)))

  if (!is.null(result$error)) {
    message("ERROR: ", result$error)
    write_csv(tibble(
      scenario = scenario, purity = purity, experiment = exp_name, note = note,
      concentration = conc, min_dist = md,
      n_clusters = NA_integer_, correct_topology = NA,
      topology = NA_character_, error = result$error
    ), file.path(out_dir, "summary_row.csv"))
    quit(status = 0)
  }

  cluster_result <- result[[1]][[1]]
  sv2cluster     <- result[[1]][[2]]
  clones         <- result[[1]][[3]]
  best_tree      <- if (!is.null(result[[2]])) result[[2]][[1]] else NULL

  n_clust  <- length(unique(sv2cluster$cluster_num))
  correct  <- is_correct_topology(best_tree)
  topo_str <- topology_string(best_tree)

  message(sprintf("  n_clusters=%d  correct=%s  topo: %s", n_clust, correct, topo_str))

  write_csv(cluster_result, file.path(out_dir, "cluster_result.csv"))
  write_csv(sv2cluster, file.path(out_dir, "sv2cluster.csv"))
  write_csv(clones,         file.path(out_dir, "cluster_centroids.csv"))
  if (!is.null(best_tree)) {
    write_csv(best_tree, file.path(out_dir, "tree_edges.csv"))
    writeLines(topo_str, file.path(out_dir, "topology.txt"))
  }

  ggsave(file.path(out_dir, "clustering_scatter.pdf"),
         make_scatter(cluster_result, scenario, purity, exp_name, conc, md),
         width = 6, height = 5.5)

  ggsave(file.path(out_dir, "clustering_scatter_delta.pdf"),
         make_delta_scatter(cluster_result, scenario, purity, exp_name, conc, md),
         width = 6, height = 5.5)

  write_csv(tibble(
    scenario = scenario, purity = purity, experiment = exp_name, note = note,
    concentration = conc, min_dist = md,
    n_clusters = n_clust, correct_topology = correct,
    topology = topo_str, error = NA_character_
  ), file.path(out_dir, "summary_row.csv"))

  message("Done.")
}

# ── MERGE ─────────────────────────────────────────────────────────────────────
run_merge <- function() {
  row_files <- list.files(SWEEP_DIR, pattern = "summary_row\\.csv",
                          recursive = TRUE, full.names = TRUE)
  if (length(row_files) == 0) stop("No summary_row.csv files found under ", SWEEP_DIR)

  message(sprintf("Merging %d summary_row.csv files ...", length(row_files)))

  summary_df <- bind_rows(lapply(row_files, read_csv, show_col_types = FALSE)) %>%
    arrange(scenario, purity, experiment, concentration, min_dist)

  out_path <- file.path(SWEEP_DIR, "summary.csv")
  write_csv(summary_df, out_path)
  message(sprintf("Written: %s  (%d rows)", out_path, nrow(summary_df)))
  message(sprintf("Correct topology: %d / %d",
                  sum(summary_df$correct_topology, na.rm = TRUE),
                  sum(!is.na(summary_df$correct_topology))))

  overview <- summary_df %>%
    group_by(scenario, purity, experiment) %>%
    summarise(
      n_correct = sum(correct_topology, na.rm = TRUE),
      n_total   = sum(!is.na(correct_topology)),
      best_conc = concentration[which.max(correct_topology)][1],
      best_md   = min_dist[which.max(correct_topology)][1],
      .groups   = "drop"
    )
  message("\nPer-case summary:")
  print(overview, n = Inf)
}

# ── DISPATCH ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || args[1] == "--help") {
  cat("Usage:\n")
  cat("  Rscript param_sweep.R --task <array_id>   # run one SLURM array task\n")
  cat("  Rscript param_sweep.R --merge             # collect results into summary.csv\n")
  quit(status = 0)
} else if (args[1] == "--task") {
  if (length(args) < 2) stop("--task requires an array ID")
  run_task(as.integer(args[2]))
} else if (args[1] == "--merge") {
  run_merge()
} else {
  stop("Unknown mode: ", args[1], "\nUse --task <id> or --merge")
}
