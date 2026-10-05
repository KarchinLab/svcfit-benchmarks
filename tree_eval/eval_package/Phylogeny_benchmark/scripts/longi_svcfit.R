#!/usr/bin/env Rscript
###############################################################################
# run_longi_pipeline.R
#
# Full longitudinal pipeline using SVCFit::run_svcfit() as the single entry
# point for all three stages: SVCFit inference, DP-GMM clustering, and tree.
#
# Stage layout
#   svcfit  : run_svcfit() for t1 and t2 independently (no clustering).
#             BEDs are saved to <svcfit_dir>/COMBAT/SVCFit_output/ so that
#             cluster_data() can locate them via data_dir.
#   cluster + tree : run_svcfit() once more with run_clustering=TRUE and
#             run_tree=TRUE (reuses t2 VCF inputs; svcfit re-run on t2 is
#             the minor cost of going through the single master function).
#
# Usage:
#   Rscript run_longi_pipeline.R \
#     --scenario S1 --purity 60 --exp exp1 \
#     --input_dir "$TREE_EVAL_LONGITUDINAL" \
#     --output_dir "$TREE_EVAL_RUN_ROOT/output" \
#     --svcfit_repo "$SVCFIT_PKG_DIR" \
#     --python_env py3 \
#     [--cov 50] [--thresh 0.1] [--stages svcfit,cluster,tree]
###############################################################################

.args <- commandArgs(trailingOnly = FALSE)
.script_file <- sub("--file=", "", .args[grep("--file=", .args)])
.script_dir  <- if (length(.script_file) > 0) dirname(normalizePath(.script_file)) else getwd()
.pub_dir     <- dirname(.script_dir)

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(reticulate)
  library(SVCFit)
  library(RColorBrewer)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(GenomicRanges)
  library(purrr)
})

option_list <- list(
  make_option(c("--scenario"), type = "character", help = "Scenario name (S1-S4)"),
  make_option(c("--purity"),   type = "integer",   help = "Tumor purity percent (e.g. 60)"),
  make_option(c("--exp"),      type = "character",  help = "Experiment name (exp1-exp5)"),
  make_option(c("--input_dir"), type = "character", help = "Existing longitudinal simulation/calling input root"),
  make_option(c("--output_dir"), type = "character", help = "New, isolated result root"),
  make_option(c("--svcfit_repo"), type = "character", help = "Pinned SVCFit Git checkout"),
  make_option(c("--python_env"), type = "character", default = "py3",
              help = "Conda environment for reticulate; set empty to use RETICULATE_PYTHON/current Python [default: py3]"),
  make_option(c("--cov"),      type = "integer",   default = 50,
              help = "Coverage [default: 50]"),
  make_option(c("--thresh"),   type = "numeric",   default = 0.1,
              help = "SVCFit CNV-order decision threshold [default: 0.1]"),
  make_option(c("--stages"),   type = "character", default = "svcfit,cluster,tree",
              help = "Comma-separated stages to run: svcfit,cluster,tree [default: all]"),
  make_option(c("--boot"),      type = "integer",   default = 0L,
              help = "Bootstrap run index (0 = original data, 1-5 = new replicates) [default: 0]"),
  make_option(c("--ccf_floor"),     type = "numeric",   default = 0.1,
              help = "Zero out CCF values below this threshold before clustering [default: 0.1]"),
  make_option(c("--concentration"), type = "numeric",   default = 1,
              help = "DP-GMM concentration parameter (higher = more clusters) [default: 1]")
)

parser <- OptionParser(option_list = option_list, add_help_option = TRUE)
opts   <- parse_args(parser)

required_opts <- c("scenario", "purity", "exp", "input_dir", "output_dir", "svcfit_repo")
missing_opts <- required_opts[vapply(required_opts, function(x) {
  value <- opts[[x]]
  is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value)
}, logical(1))]
if (length(missing_opts)) stop("Missing required option(s): ", paste(missing_opts, collapse = ", "))
if (!dir.exists(opts$input_dir)) stop("Input directory does not exist: ", opts$input_dir)
if (!dir.exists(opts$svcfit_repo)) stop("SVCFit checkout does not exist: ", opts$svcfit_repo)

# Select Python explicitly when requested. RETICULATE_PYTHON takes precedence;
# an empty --python_env leaves selection to the active environment/reticulate.
reticulate_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
if (nzchar(reticulate_python)) {
  use_python(reticulate_python, required = TRUE)
} else if (nzchar(opts$python_env)) {
  use_condaenv(opts$python_env, required = TRUE)
}

stages <- trimws(strsplit(opts$stages, ",")[[1]])

atomic_write <- function(path, writer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  writer(tmp)
  if (!file.exists(tmp) || file.info(tmp)$size <= 0)
    stop("Atomic writer produced an empty file: ", path)
  if (!file.rename(tmp, path)) stop("Could not atomically install output: ", path)
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  atomic_write(path, function(tmp) write_csv(x, tmp))
}

atomic_write_delim <- function(x, path) {
  atomic_write(path, function(tmp) write_delim(x, tmp, delim = "\t", quote = "none"))
}

atomic_write_lines <- function(x, path) {
  atomic_write(path, function(tmp) writeLines(x, tmp))
}

atomic_save_rds <- function(x, path) {
  atomic_write(path, function(tmp) saveRDS(x, tmp))
}

# BOOT=0: original path structure (backward compatible)
# BOOT>0: b{boot}/ subdirectory for new bootstrap replicates
if (opts$boot == 0L) {
  source_longi_dir <- file.path(opts$input_dir, opts$scenario,
                                paste0("c", opts$cov, "p", opts$purity), opts$exp)
  longi_dir <- file.path(opts$output_dir, opts$scenario,
                         paste0("c", opts$cov, "p", opts$purity), opts$exp)
} else {
  source_longi_dir <- file.path(opts$input_dir, opts$scenario,
                                paste0("c", opts$cov, "p", opts$purity),
                                paste0("b", opts$boot), opts$exp)
  longi_dir <- file.path(opts$output_dir, opts$scenario,
                         paste0("c", opts$cov, "p", opts$purity),
                         paste0("b", opts$boot), opts$exp)
}
if (!dir.exists(source_longi_dir)) stop("Case input directory does not exist: ", source_longi_dir)
dir.create(longi_dir, recursive = TRUE, showWarnings = FALSE)
svcfit_dir <- file.path(longi_dir, "svcfit_output")
# cluster_data() looks for <data_dir>/COMBAT/SVCFit_output/<sample>.bed
combat_dir <- file.path(svcfit_dir, "COMBAT", "SVCFit_output")

# True per-timepoint purity from simulation design (mirrors longi_short.sh).
# Use ground-truth purity, NOT FACETS-estimated purity.
# t1 purity = ppur for all scenarios; t2 purity differs only in S4 (purity halves).
pur_frac <- opts$purity / 100
true_purity <- list(
  S1 = c(t1 = pur_frac,           t2 = pur_frac),          # subclone swap:       purity constant
  S2 = c(t1 = pur_frac,           t2 = pur_frac),          # negative control:    purity constant
  S3 = c(t1 = pur_frac,           t2 = pur_frac),          # gradual drift:       purity constant
  S4 = c(t1 = pur_frac,           t2 = pur_frac / 2)       # asymmetric regression: purity halves at t2
)
purity_t1 <- unname(true_purity[[opts$scenario]]["t1"])
purity_t2 <- unname(true_purity[[opts$scenario]]["t2"])

# Helper: build per-timepoint input paths
tp_paths <- function(tp) {
  list(
    het_snp = file.path(source_longi_dir, tp, "SNP",
                        sprintf("het_near_sv_%s_%s.vcf", opts$exp, tp)),
    on_sv   = file.path(source_longi_dir, tp, "SNP",
                        sprintf("het_on_sv_%s_%s.vcf",  opts$exp, tp)),
    sv      = file.path(source_longi_dir, tp, "svtyp",
                        sprintf("svt_%s_%s.vcf",         opts$exp, tp)),
    cnv     = file.path(source_longi_dir, tp, "facet",
                        sprintf("%s_%s.bed",              opts$exp, tp))
  )
}

###############################################################################
# Stage 1: SVCFit inference for t1 and t2
###############################################################################
if ("svcfit" %in% stages) {
  dir.create(combat_dir, recursive = TRUE, showWarnings = FALSE)

  for (tp in c("t1", "t2")) {
    message(sprintf("\n=== SVCFit: %s / %s / %s ===", opts$scenario, opts$exp, tp))
    p <- tp_paths(tp)

    samp_name <- sprintf("%s_p%s_%s_%s", opts$scenario, opts$purity, opts$exp, tp)

    result <- run_svcfit(
      p_het  = p$het_snp,
      p_onsv = p$on_sv,
      p_sv   = p$sv,
      p_cnv  = p$cnv,
      samp   = samp_name,
      exper  = opts$exp,
      thresh = opts$thresh
    )

    # Save to COMBAT layout (required by cluster_data) and to svcfit_dir directly
    bed_combat <- file.path(combat_dir,  sprintf("%s.bed", samp_name))
    bed_flat   <- file.path(svcfit_dir,  sprintf("%s.bed", samp_name))
    atomic_write_delim(result$svcf, bed_combat)
    atomic_write_delim(result$svcf, bed_flat)
    message(sprintf("Written %d SVs to %s", nrow(result$svcf), bed_combat))
  }
}

###############################################################################
# Stage 2+3: Clustering and tree via run_svcfit()
#
# run_svcfit() with run_clustering=TRUE reads both BEDs via cluster_data()
# using data_dir = svcfit_dir (→ svcfit_dir/COMBAT/SVCFit_output/<sample>.bed).
# run_tree=TRUE then calls build_tree() on the cluster CCF table.
###############################################################################
if (any(c("cluster", "tree") %in% stages)) {

  run_tree <- "tree" %in% stages

  # A tree result is valid only for the clustering created by this invocation.
  # Clear previous downstream outcomes before clustering starts, then write
  # exactly one of tree_result.rds or NO_TREE.txt below.
  cluster_dir <- file.path(longi_dir, "clustering")
  tree_dir <- file.path(longi_dir, "tree")
  if (dir.exists(cluster_dir)) unlink(cluster_dir, recursive = TRUE, force = TRUE)
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  if (run_tree && dir.exists(tree_dir)) unlink(tree_dir, recursive = TRUE, force = TRUE)
  if (run_tree) dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)

  # Write pair_path: no-header TSV with pre_BAT and on_BAT columns
  pair_path <- file.path(svcfit_dir, "pair_path.txt")
  write.table(
    data.frame(pre_BAT = sprintf("%s_p%s_%s_t1", opts$scenario, opts$purity, opts$exp),
               on_BAT  = sprintf("%s_p%s_%s_t2", opts$scenario, opts$purity, opts$exp)),
    pair_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  # Write pur_path: TSV with 'sample' and 'purity' columns
  pur_path <- file.path(svcfit_dir, "pur_path.txt")
  write.table(
    data.frame(sample  = sprintf("%s_p%s_%s_%s", opts$scenario, opts$purity, opts$exp, c("t1", "t2")),
               purity  = c(purity_t1, purity_t2)),
    pur_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  message(sprintf("\n=== Clustering%s: %s / %s ===",
                  if (run_tree) " + tree" else "", opts$scenario, opts$exp))

  full_result <- build_trees(
    run_clustering            = TRUE,
    pair_path                 = pair_path,
    pur_path                  = pur_path,
    data_dir                  = svcfit_dir,
    pair_num                  = 1L,
    deduplicate               = TRUE,
    concentration             = 1,
    min_dist                  = 0.2,

    run_tree                  = run_tree,
    lineage_precedence_thresh = 0.2,
    sum_filter_thresh         = 0.2,
    ccf_floor                 = opts$ccf_floor,
    linear_penalty            = 0.3
  )

  # Save clustering outputs
  cluster_result <- full_result[[1]][[1]]   # full cluster-assignment table
  sv2cluster     <- full_result[[1]][[2]]   # SV-to-cluster mapping
  clones         <- full_result[[1]][[3]]   # clone CCF table for pair 1

  atomic_write_csv(cluster_result, file.path(cluster_dir, "cluster_result.csv"))
  atomic_write_csv(sv2cluster,     file.path(cluster_dir, "sv2cluster.csv"))
  atomic_write_csv(clones,         file.path(cluster_dir, "cluster_centroids.csv"))
  message(sprintf("Clustering output saved to %s", cluster_dir))

  # Save exactly one explicit tree outcome.
  if (run_tree && !is.null(full_result[[2]])) {
    best_tree <- full_result[[2]][[1]]
    mcf_mat   <- full_result[[2]][[2]]

    topo_str <- paste(apply(best_tree, 1,
                            function(r) paste(r["parent"], r["child"], sep = "->")),
                      collapse = "; ")

    atomic_write_csv(best_tree, file.path(tree_dir, "tree_edges.csv"))
    atomic_write_lines(topo_str, file.path(tree_dir, "topology.txt"))
    atomic_save_rds(full_result[[2]], file.path(tree_dir, "tree_result.rds"))
    message(sprintf("Tree output saved to %s", tree_dir))
    message(sprintf("Topology: %s", topo_str))
  } else if (run_tree) {
    no_tree <- file.path(tree_dir, "NO_TREE.txt")
    atomic_write_lines(c(
      "NO_TREE",
      sprintf("scenario=%s", opts$scenario),
      sprintf("purity=%s", opts$purity),
      sprintf("boot=%s", opts$boot),
      sprintf("experiment=%s", opts$exp),
      "reason=build_trees returned no admissible topology"
    ), no_tree)
    message(sprintf("No admissible topology; marker written to %s", no_tree))
  }
}

message(sprintf("\n=== Pipeline complete: %s / %s ===", opts$scenario, opts$exp))
