#!/usr/bin/env Rscript
#
# A2 coverage correlation — per-run AO extraction.
#
# Usage:
#   Rscript extract_per_run.R <purity> <replicate> <timepoint>
#
# For one (purity, replicate, timepoint):
#   1. Parses the SVtyper VCF and extracts AO from the `bulk` sample column
#      for each SV.
#   2. Reads the inferred-cluster assignment CSV (sv_id -> inferred_cluster_id).
#   3. Computes mean AO per inferred cluster.
#   4. Writes one CSV with one row per inferred cluster, including metadata.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript extract_per_run.R <purity> <replicate> <timepoint>")
}
purity    <- as.integer(args[1])
replicate <- as.integer(args[2])
timepoint <- args[3]
stopifnot(timepoint %in% c("pre", "post"))

script_dir <- Sys.getenv("A2_SCRIPT_DIR", unset = ".")
source(file.path(script_dir, "config.R"))

dir.create(PER_RUN_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Parse VCF, extract AO from `bulk` column.
# -----------------------------------------------------------------------------

vcf_file <- vcf_path(purity, replicate, timepoint)
if (!file.exists(vcf_file)) stop(sprintf("VCF not found: %s", vcf_file))

vcf_lines  <- readLines(vcf_file, warn = FALSE)
header_idx <- grep("^#CHROM", vcf_lines)
if (length(header_idx) != 1L) stop("Expected exactly one #CHROM header line.")
header_cols <- strsplit(sub("^#", "", vcf_lines[header_idx]), "\t", fixed = TRUE)[[1]]
data_lines  <- vcf_lines[!startsWith(vcf_lines, "#")]
data_lines  <- data_lines[nzchar(data_lines)]

# The example VCF has two identical `bulk` columns; take the first.
bulk_idx   <- which(header_cols == "bulk")
if (length(bulk_idx) == 0L) stop("No `bulk` sample column in VCF header.")
bulk_idx   <- bulk_idx[1]
format_idx <- which(header_cols == "FORMAT")
id_idx     <- which(header_cols == "ID")

extract_one <- function(line) {
  fields  <- strsplit(line, "\t", fixed = TRUE)[[1]]
  fmt_keys <- strsplit(fields[format_idx], ":", fixed = TRUE)[[1]]
  ao_idx   <- which(fmt_keys == "AO")
  if (length(ao_idx) == 0L) {
    return(c(sv_id = fields[id_idx], AO = NA_character_))
  }
  bulk_vals <- strsplit(fields[bulk_idx], ":", fixed = TRUE)[[1]]
  c(sv_id = fields[id_idx], AO = bulk_vals[ao_idx])
}

vcf_ao <- as.data.frame(do.call(rbind, lapply(data_lines, extract_one)),
                        stringsAsFactors = FALSE)
vcf_ao$AO <- suppressWarnings(as.numeric(vcf_ao$AO))

n_no_ao <- sum(is.na(vcf_ao$AO))
if (n_no_ao > 0L) {
  message(sprintf("Note: %d / %d VCF records have no parseable AO.",
                  n_no_ao, nrow(vcf_ao)))
}

# -----------------------------------------------------------------------------
# 2. Read inferred-cluster assignments (one file covers both timepoints).
# -----------------------------------------------------------------------------

ca_file <- cluster_assign_path(purity, replicate)
if (!file.exists(ca_file)) stop(sprintf("Cluster assignments not found: %s", ca_file))
clust_full <- read.csv(ca_file, stringsAsFactors = FALSE)
required <- c("shared", "cluster_num", "stage")
missing  <- setdiff(required, colnames(clust_full))
if (length(missing) > 0L) {
  stop(sprintf("Cluster file %s missing columns: %s",
               ca_file, paste(missing, collapse = ", ")))
}

# Filter to the requested timepoint.
# For shared SVs, `shared` holds the VCF ID.  For timepoint-unique SVs (shared
# is NA), the VCF ID is embedded at the end of `event_id` (e.g.
# "pair1_S1_p20_exp1_t1_MantaDEL:5:0:1:0:0:0" -> "MantaDEL:5:0:1:0:0:0").
stage_val <- if (timepoint == "pre") "pre_BAT" else "on_BAT"
clust <- clust_full[clust_full$stage == stage_val, ]
sv_id <- ifelse(
  !is.na(clust$shared),
  clust$shared,
  sub("^.*_(Manta.+)$", "\\1", clust$event_id)
)
clust <- data.frame(sv_id = sv_id, inferred_cluster_id = clust$cluster_num,
                    stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 3. Join AO into cluster assignment, compute per-cluster mean AO.
# -----------------------------------------------------------------------------

joined <- merge(clust, vcf_ao, by = "sv_id", all.x = TRUE)
n_unjoined <- sum(is.na(joined$AO))
if (n_unjoined > 0L) {
  message(sprintf("Note: %d / %d clustered SVs had no AO match in the VCF (likely BND mates dropped by A2-fix1).",
                  n_unjoined, nrow(joined)))
}

per_cluster <- aggregate(
  AO ~ inferred_cluster_id,
  data = joined,
  FUN  = function(x) c(n_svs = length(x), mean_AO = mean(x, na.rm = TRUE))
)
per_cluster <- do.call(data.frame, per_cluster)
colnames(per_cluster) <- c("inferred_cluster_id", "n_svs", "mean_AO")

# -----------------------------------------------------------------------------
# 4. Add metadata, write.
# -----------------------------------------------------------------------------

per_cluster$purity    <- purity
per_cluster$replicate <- replicate
per_cluster$timepoint <- timepoint
per_cluster <- per_cluster[, c("purity", "replicate", "timepoint",
                               "inferred_cluster_id", "n_svs", "mean_AO")]

out_path <- file.path(
  PER_RUN_DIR,
  sprintf("eff_cov_p%02d_r%02d_%s.csv", purity, replicate, timepoint)
)
write.csv(per_cluster, out_path, row.names = FALSE)
message(sprintf("Wrote %s  (%d clusters)", out_path, nrow(per_cluster)))
