#!/usr/bin/env bash
# Cluster-neutral configuration template.
#
# Copy this file to config.local.sh and set PROJECT_ROOT plus
# SVCFIT_RUNTIME_LOADER for the local cluster. Site-specific scheduler options
# are optional; workflow entrypoints retain their existing names and arguments.

: "${PROJECT_ROOT:?Set PROJECT_ROOT to the numbered project root}"
: "${SVCFIT_RUNTIME_LOADER:?Set SVCFIT_RUNTIME_LOADER to the site runtime loader}"

export SOFTWARE_ROOT="${SOFTWARE_ROOT:-$PROJECT_ROOT/01_software}"
export DATA_ROOT="${DATA_ROOT:-$PROJECT_ROOT/02_data}"
export ANALYSIS_ROOT="${ANALYSIS_ROOT:-$PROJECT_ROOT/03_analysis}"
export QC_ROOT="${QC_ROOT:-$PROJECT_ROOT/04_qc}"
export PUBLICATION_ROOT="${PUBLICATION_ROOT:-$PROJECT_ROOT/05_publication}"

[[ -r "$SVCFIT_RUNTIME_LOADER" ]] || {
    echo "ERROR: missing site runtime loader: $SVCFIT_RUNTIME_LOADER" >&2
    return 1 2>/dev/null || exit 1
}
# shellcheck disable=SC1090
source "$SVCFIT_RUNTIME_LOADER"

export VISOR_ROOT="${VISOR_ROOT:-$SOFTWARE_ROOT/svcfit-benchmarks}"
export CHRX_DIR="${CHRX_DIR:-$ANALYSIS_ROOT/visor_chrx/runs}"
export REPLICATES_DIR="${REPLICATES_DIR:-$ANALYSIS_ROOT/visor/runs}"
export TREE_EVAL_DIR="${TREE_EVAL_DIR:-$VISOR_ROOT/tree_eval}"
export TREE_EVAL_LONGITUDINAL="${TREE_EVAL_LONGITUDINAL:-$DATA_ROOT/tree_eval/longitudinal}"
export TREE_EVAL_TRUTH_DIR="${TREE_EVAL_TRUTH_DIR:-$DATA_ROOT/tree_eval/input_data/hack}"
export REF_CHRX="${REF_CHRX:-$DATA_ROOT/references/genomes/chr22-X.fa}"
export REF_AUTO="${REF_AUTO:-$DATA_ROOT/tree_eval/input_data/reference/chr1-2.fa}"

export SVCFIT_ROOT="${SVCFIT_ROOT:-$PROJECT_ROOT}"
export MC_HACK_BASE="${MC_HACK_BASE:-$DATA_ROOT/visor/truth}"
export NORM_SHORT_DIR="${NORM_SHORT_DIR:-$DATA_ROOT/visor/norm_short}"
export SVCFIT_RESOURCES="${SVCFIT_RESOURCES:-$DATA_ROOT/references/svcfit}"
export SVCFIT_PKG_DIR="${SVCFIT_PKG_DIR:-$SOFTWARE_ROOT/SVCFit}"
export SVCFIT_SRC="${SVCFIT_SRC:-$SVCFIT_PKG_DIR}"
export SVCFIT_PATH="${SVCFIT_PATH:-$SVCFIT_PKG_DIR}"

export PROSTATE_DATA_DIR="${PROSTATE_DATA_DIR:-$DATA_ROOT/prostate_mixture}"
export PROSTATE_REF="${PROSTATE_REF:-$PROSTATE_DATA_DIR/ref/hs37d5.fa}"
export PIPELINE_CONDA_SH="${PIPELINE_CONDA_SH:-${CONDA_SH:-}}"
export VISOR_ENV="${VISOR_ENV:-${ENV_VISOR:-}}"
export MANTA_HOME="${MANTA_HOME:-${ENV_MANTA:-}}"
export MANTA_CONFIG="${MANTA_CONFIG:-$MANTA_HOME/bin/configManta.py}"
if [[ -z "${MANTA_CONVERT:-}" ]]; then
    if [[ -x "$MANTA_HOME/bin/convertInversion.py" ]]; then
        export MANTA_CONVERT="$MANTA_HOME/bin/convertInversion.py"
    else
        export MANTA_CONVERT="$MANTA_HOME/libexec/convertInversion.py"
    fi
fi
export FACET_HOME="${FACET_HOME:-$SOFTWARE_ROOT/facets}"
export FACET_SNP_PILEUP="${FACET_SNP_PILEUP:-$FACET_HOME/snp-pileup}"
export FACET_R_SCRIPT="${FACET_R_SCRIPT:-$FACET_HOME/facet.R}"
export SAMTOOLS="${SAMTOOLS:-${ENV_VISOR:+$ENV_VISOR/bin/samtools}}"
# R with DNAcopy for the chromosome X depth segmentation (visor_chrX 06, prostate 04a);
# built from tools/environments/dnacopy.yml by tools/setup_rockfish_tool_envs.sh.
export ENV_DNACOPY="${ENV_DNACOPY:-}"
export CHRX_RSCRIPT="${CHRX_RSCRIPT:-${ENV_DNACOPY:+$ENV_DNACOPY/bin/Rscript}}"

export SHORT_DIR="${SHORT_DIR:-$CHRX_DIR/short}"
export CALLS_DIR="${CALLS_DIR:-$CHRX_DIR/calls}"
export ENV_SVCFIT="${ENV_SVCFIT:-}"
export RLIB="${RLIB:-$ENV_SVCFIT/lib/R/library}"
export NPJ_PYTHON_ENV="${NPJ_PYTHON_ENV:-$ENV_SVCFIT}"
export SVCFIT_R="${SVCFIT_R:-${RSCRIPT:-$ENV_SVCFIT/bin/Rscript}}"
export SVCFIT_RUNTIME_LOCKFILE="${SVCFIT_RUNTIME_LOCKFILE:-$ENV_SVCFIT/conda-explicit-linux-64.txt}"

# Optional site scheduler settings. Leave unset to use cluster defaults.
# Only run_svcfit_and_evaluate.sh and tools/run_svcfit_self_evaluation.sh read
# all three (25_chrx_rescore_svcfit.sh reads SLURM_PARTITION). For the other
# launchers, export sbatch's own SBATCH_ACCOUNT, SBATCH_PARTITION and
# SBATCH_QOS, or edit the #SBATCH headers.
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"
export SLURM_PARTITION="${SLURM_PARTITION:-}"
export SLURM_QOS="${SLURM_QOS:-}"
export CHRX_EXCLUDE_NODES="${CHRX_EXCLUDE_NODES:-}"
if [[ -n "${SLURM_BIN:-}" ]]; then
    export PATH="$SLURM_BIN:$PATH"
fi

_svcfit_require() {
    local name val
    for name in "$@"; do
        val="${!name-}"
        [[ -n "$val" ]] || { echo "config.local.sh: $name is unset" >&2; return 1; }
        [[ -e "$val" ]] || { echo "config.local.sh: $name -> $val does not exist" >&2; return 1; }
    done
}
_svcfit_require PROJECT_ROOT SOFTWARE_ROOT DATA_ROOT ANALYSIS_ROOT QC_ROOT \
    VISOR_ROOT CHRX_DIR REPLICATES_DIR TREE_EVAL_LONGITUDINAL TREE_EVAL_TRUTH_DIR \
    REF_CHRX REF_AUTO SVCFIT_PKG_DIR PROSTATE_DATA_DIR PROSTATE_REF \
    CONDA_SH RETICULATE_PYTHON ENV_SVCFIT ENV_VISOR ENV_GATK ENV_FACET \
    ENV_MANTA ENV_SVTYPER ENV_SVCLONE MANTA_CONFIG MANTA_CONVERT FACET_HOME \
    FACET_SNP_PILEUP FACET_R_SCRIPT SHORT_DIR CALLS_DIR \
    || { echo "Fix config.local.sh before running anything." >&2; return 1 2>/dev/null || exit 1; }

export VISOR_CONFIG_LOADED=1
unset -f _svcfit_require
