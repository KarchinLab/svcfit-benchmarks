#!/usr/bin/env bash
# Bundles all assets for a self-contained, publishable experiment package.
# BAM and BAI files are intentionally excluded.
#
# After running this script, the directory layout will be:
#   publishable/
#     data/
#       reference/     hg38 chr1-2 FASTA + BWA/samtools indices
#       hack/          VISOR HACk haplotype FASTAs + clone SV BED files (e1-e5, snp1, snp1_1)
#       norm_short/    Normal sample BED descriptors (no BAMs)
#       SVCFit/        SVCFit R package source
#     scripts/         All pipeline and analysis scripts (+ a2_coverage_correlation/)
#     figures/         Produced plots and figures (.pdf, .png, .svg)
#     outputs/         SVCFit results and downstream analysis outputs
#     design/          Experiment design documents, READMEs, and protocol notes
#
# Paths in the pipeline scripts that will need updating after the copy:
#   longi_short.sh   : hack_base, ref
#   longi_calling.sh : ref, nbam (norm BAM path; BAMs are NOT copied — supply your own)
#   longi_svcfit.R   : devtools::load_all() path for SVCFit, PYTHON_ENV
#   a2_coverage_correlation/config.R : any absolute paths
#
# Usage: bash centralize_experiment.sh [DEST_DIR]
#   DEST_DIR defaults to this script's directory (VISOR_longitudinal)

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${1:-$SCRIPT_DIR}"

# --- machine-specific paths ---------------------------------------------------
# Paths come from config.local.sh, found by walking UP from this file -- the same
# mechanism visor_chrX/scripts/chrx_common.sh uses. $VISOR_CONFIG overrides it.
# See README.md.
if [[ -z "${VISOR_CONFIG_LOADED:-}" ]]; then
    _cfg="${VISOR_CONFIG:-}"
    if [[ -z "$_cfg" ]]; then
        _d="$SCRIPT_DIR"
        while [[ "$_d" != "/" ]]; do
            [[ -r "$_d/config.local.sh" ]] && { _cfg="$_d/config.local.sh"; break; }
            _d="$(dirname "$_d")"
        done
    fi
    [[ -n "$_cfg" && -r "$_cfg" ]] || {
        echo "ERROR: config.local.sh not found by walking up from ${BASH_SOURCE[0]}." >&2
        echo "       cp config.example.sh config.local.sh at the repo root and edit it." >&2
        echo "       See README.md." >&2
        exit 1
    }
    # shellcheck disable=SC1090
    source "$_cfg"
    unset _cfg _d
fi

# --- SOURCE PATHS (remote data) ---
# REF_AUTO is the combined chr1-2 FASTA; the package wants the per-chromosome
# directory that holds it, which is its parent.
SRC_REF="$(dirname "$REF_AUTO")"
SRC_HACK="$MC_HACK_BASE"
SRC_NORM="$NORM_SHORT_DIR"
SRC_SVCFIT="$SVCFIT_PKG_DIR"

# --- SOURCE PATHS (local — outputs and figures produced by the pipeline) ---
# SVCFit pipeline outputs used by evaluate_downstream.R (--work_dir argument).
# Contains S1/c50p{purity}/[b{boot}/]{exp}/{clustering,tree,t1,t2}/
#
# This is the 5.9 TB raw simulation tree, deliberately NOT moved into the
# repository when the rest of tree_eval was, so it keeps its own variable rather
# than being derived from this script's location. See README.md.
: "${TREE_EVAL_LONGITUDINAL:?set TREE_EVAL_LONGITUDINAL in config.local.sh -- the raw longitudinal simulation tree}"
SRC_OUTPUTS="$TREE_EVAL_LONGITUDINAL"

# --- DESTINATION DIRS ---
DEST_PUB="$DEST/publishable"
DEST_DATA="$DEST_PUB/data"
DEST_REF="$DEST_DATA/reference"
DEST_HACK="$DEST_DATA/hack"
DEST_NORM="$DEST_DATA/norm_short"
DEST_SVCFIT="$DEST_DATA/SVCFit"
DEST_SCRIPTS="$DEST_PUB/scripts"
DEST_FIGURES="$DEST_PUB/figures"
DEST_OUTPUTS="$DEST_PUB/outputs"
DEST_DESIGN="$DEST_PUB/design"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
hr()  { echo "────────────────────────────────────────────────────────────"; }

hr
log "Building publishable package into: $DEST_PUB"
hr

# Verify remote data sources exist
missing=0
for src in "$SRC_REF" "$SRC_HACK" "$SRC_NORM" "$SRC_SVCFIT"; do
    if [[ ! -d "$src" ]]; then
        log "ERROR: Source not found: $src"
        missing=1
    fi
done
(( missing )) && { log "Aborting — fix missing sources above."; exit 1; }

mkdir -p "$DEST_REF" "$DEST_HACK" "$DEST_NORM" "$DEST_SVCFIT" \
         "$DEST_SCRIPTS" "$DEST_FIGURES" "$DEST_OUTPUTS" "$DEST_DESIGN"

# ── 1. Reference genome ──────────────────────────────────────────────────────
log "Copying reference genome (chr1-2.fa + indices) …"
rsync -av --progress \
    --include="chr1-2.fa" \
    --include="chr1-2.fa.*" \
    --include="chr1-2.dict" \
    --exclude="*" \
    "$SRC_REF/" "$DEST_REF/"
log "Reference: done."
hr

# ── 2. VISOR HACk haplotype genomes (no BAMs) ────────────────────────────────
log "Copying VISOR HACk genomes (e1–e5, snp1, snp1_1) + clone BEDs …"
log "  Source size: $(du -sh "$SRC_HACK" 2>/dev/null | cut -f1)"
rsync -av --progress \
    --exclude="*.bam" \
    --exclude="*.bai" \
    "$SRC_HACK/" "$DEST_HACK/"
log "HACk genomes: done."
hr

# ── 3. Normal sample BED descriptors (no BAMs) ───────────────────────────────
log "Copying normal sample BED descriptors (BAMs excluded) …"
rsync -av --progress \
    --exclude="*.bam" \
    --exclude="*.bai" \
    "$SRC_NORM/" "$DEST_NORM/"
log "norm_short: done."
hr

# ── 4. SVCFit R package source ───────────────────────────────────────────────
log "Copying SVCFit R package source …"
rsync -av --progress \
    --exclude=".git" \
    "$SRC_SVCFIT/" "$DEST_SVCFIT/"
log "SVCFit: done."
hr

# ── 5. Scripts ────────────────────────────────────────────────────────────────
log "Copying pipeline and analysis scripts …"
rsync -av --progress \
    --include="*.sh" \
    --include="*.R" \
    --include="*.py" \
    --exclude="*" \
    "$SCRIPT_DIR/" "$DEST_SCRIPTS/"
# Copy a2_coverage_correlation subdirectory (scripts + configs)
rsync -av --progress \
    --exclude="*.bam" \
    --exclude="*.bai" \
    "$SCRIPT_DIR/a2_coverage_correlation/" "$DEST_SCRIPTS/a2_coverage_correlation/"
log "Scripts: done."
hr

# ── 6. Figures ────────────────────────────────────────────────────────────────
log "Copying figures (.pdf, .png, .svg) …"
rsync -av --progress \
    --include="*.pdf" \
    --include="*.png" \
    --include="*.svg" \
    --exclude="*" \
    "$SCRIPT_DIR/" "$DEST_FIGURES/"
log "Figures: done."
hr

# ── 7. Outputs (S1 pipeline results needed by evaluate_downstream.R) ──────────
if [[ -n "$SRC_OUTPUTS" && -d "$SRC_OUTPUTS" ]]; then
    log "Copying S1 pipeline outputs from $SRC_OUTPUTS …"
    log "  Source size: $(du -sh "$SRC_OUTPUTS/S1" 2>/dev/null | cut -f1)"
    # Only S1 (the scenario evaluate_downstream.R reads).
    # Per experiment keeps: clustering/ CSVs, tree/ RDS, t{1,2}/svtyp/ VCFs,
    # t{1,2}/facet/ BED files.
    # Excludes: BAMs/BAIs, large FACET intermediates (.csv.gz, .RData),
    #           SNP/, manta/, short.out/, svcfit_output/ (not read by the script).
    rsync -av --progress \
        --exclude="*.bam" \
        --exclude="*.bai" \
        --exclude="*.csv.gz" \
        --exclude="*.RData" \
        --exclude="*.pickle" \
        --exclude="short.out/" \
        --exclude="short.bed" \
        --exclude="SNP/" \
        --exclude="manta/" \
        --exclude="svcfit_output/" \
        "$SRC_OUTPUTS/S1/" "$DEST_OUTPUTS/S1/"
    log "Outputs: done."
else
    log "SRC_OUTPUTS not set or not found — skipping outputs section."
fi
hr

# ── 8. Experiment design documents ───────────────────────────────────────────
log "Copying experiment design documents …"
rsync -av --progress \
    --include="*.md" \
    --include="README" \
    --exclude="*" \
    "$SCRIPT_DIR/" "$DEST_DESIGN/"
log "Design docs: done."
hr

# ── 9. Summary ───────────────────────────────────────────────────────────────
log "All assets copied.  Final layout:"
du -sh "$DEST_PUB"/*/
hr

cat <<'EOF'

NEXT STEPS — update hardcoded paths in the scripts to point to the centralized data:

  publishable/scripts/longi_short.sh
    hack_base=<DEST>/publishable/data     (directory containing hack/)
    ref=<DEST>/publishable/data/reference/chr1-2.fa

  publishable/scripts/longi_calling.sh
    ref=<DEST>/publishable/data/reference/chr1-2.fa
    nbam= ... (BAMs were not copied — supply the path to your normal BAMs separately)

  publishable/scripts/longi_svcfit.R
    devtools::load_all('<DEST>/publishable/data/SVCFit')

  publishable/scripts/a2_coverage_correlation/config.R
    Review all absolute paths and update accordingly.

  evaluate_downstream.R
    --work_dir  <DEST>/publishable/outputs   (contains S1/)
    --truth_dir <DEST>/publishable/data/hack
    --out_dir   <DEST>/publishable/outputs/evaluation

EOF
log "Done."
