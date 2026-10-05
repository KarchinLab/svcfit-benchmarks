# A2 coverage correlation — configuration
#
# Sourced by extract_per_run.R and aggregate_and_analyze.R.
# Edit the FILL IN comments before running.

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

# Root of A2-simple simulation outputs.
SIM_ROOT <- file.path(dirname(dirname(normalizePath(Sys.getenv("A2_SCRIPT_DIR", unset = ".")))), "outputs", "S1")

# Existing signed-error table from the Fig A2-1B pipeline.
# REQUIRED columns: purity, replicate, timepoint, true_clone,
#                   inferred_cluster_id, true_ccf, inferred_ccf, signed_error
# - purity: integer (20 / 40 / 60 / 80)
# - replicate: integer (1..20)
# - timepoint: "pre" or "post"
# - true_clone: "Trunk" / "Subclone 1" / "Subclone 2"
# - inferred_cluster_id: matches IDs used in cluster_assign_path() output
# If the Fig A2-1B pipeline does not currently export this as a single CSV,
# you'll need to export it once from whatever step produces panel B.
SIGNED_ERROR_TABLE <- file.path(SIM_ROOT, "A2_simple", "signed_errors.csv")   # FILL IN

# Path to per-run SVtyper VCF (the file that contains the `bulk` sample column
# and `AO` FORMAT subfield, as in the example svt_exp1_t1.vcf).
vcf_path <- function(purity, rep, timepoint) {
  tpoint_num <- if (timepoint == "pre") 1L else 2L
  file.path(
    SIM_ROOT,
    sprintf("c50p%d", purity),
    sprintf("exp%d", rep),
    sprintf("t%d", tpoint_num),
    "svtyp",
    sprintf("svt_exp%d_t%d.vcf", rep, tpoint_num)
  )
}

# Path to the per-run inferred-cluster assignment CSV.
# This file covers BOTH timepoints; filter by `stage` in extract_per_run.R.
# REQUIRED columns: shared (VCF ID), cluster_num, stage (pre_BAT / on_BAT).
cluster_assign_path <- function(purity, rep) {
  file.path(
    SIM_ROOT,
    sprintf("c50p%d", purity),
    sprintf("exp%d", rep),
    "clustering", "sv2cluster.csv"
  )
}

# -----------------------------------------------------------------------------
# Sweep grid (must match A2-simple)
# -----------------------------------------------------------------------------

PURITIES   <- c(20L, 40L, 60L, 80L)
REPLICATES <- 1L:20L
TIMEPOINTS <- c("pre", "post")

# -----------------------------------------------------------------------------
# Output locations
# -----------------------------------------------------------------------------

OUT_DIR     <- file.path(SIM_ROOT, "A2_coverage_correlation")
PER_RUN_DIR <- file.path(OUT_DIR, "per_run")

# -----------------------------------------------------------------------------
# Bootstrap settings
# -----------------------------------------------------------------------------

BOOTSTRAP_B    <- 1000L
BOOTSTRAP_SEED <- 42L
