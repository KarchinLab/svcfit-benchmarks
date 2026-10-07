#!/bin/bash
set -euo pipefail
# Shared variables for the VISOR replicate pipeline.
# 30 replicates x 5 experiments x 15 conditions = 2250 array tasks (0-2249)
#
# Encoding:
#   rep_idx  = SLURM_ARRAY_TASK_ID / 75    (0-29 → rep1-rep30)
#   exp_idx  = (SLURM_ARRAY_TASK_ID % 75) / 15   (0-4 → exp1-exp5)
#   cond_idx = SLURM_ARRAY_TASK_ID % 15    (0-14 → 15 conditions)
#
# Condition index → (purity, mixture):
#   0: p10m10  1: p10m30  2: p10m50
#   3: p20m10  4: p20m30  5: p20m50
#   6: p40m10  7: p40m30  8: p40m50
#   9: p60m10 10: p60m30 11: p60m50
#  12: p80m10 13: p80m30 14: p80m50

# --- Decode task ID ---
rep_idx=$(( SLURM_ARRAY_TASK_ID / 75 ))
remainder=$(( SLURM_ARRAY_TASK_ID % 75 ))
exp_idx=$(( remainder / 15 ))
cond_idx=$(( remainder % 15 ))

rep=rep$(( rep_idx + 1 ))

# --- Experiment mapping ---
exp_names=(exp1 exp2 exp3 exp4 exp5)
hack_names=(e1 e2 e3 e4 e5)
snp_names=(snp1 snp1_1 snp1_1 snp1_1 snp1_1)
normal_names=(normal o_normal o_normal o_normal o_normal)

exp_name=${exp_names[$exp_idx]}
hack_name=${hack_names[$exp_idx]}
snp_name=${snp_names[$exp_idx]}
normal_name=${normal_names[$exp_idx]}

# --- Condition mapping ---
purities=(10 10 10 20 20 20 40 40 40 60 60 60 80 80 80)
mixtures=(10 30 50 10 30 50 10 30 50 10 30 50 10 30 50)

ppur=${purities[$cond_idx]}
sub_mix=${mixtures[$cond_idx]}
cov=50
pur=100

# Condition name: c50p{ppur}m{sub_mix}
cond_name=c${cov}p${ppur}m${sub_mix}

# Clone fractions: c2 (subclone), c3 (major clone), normal
# sub_mix% of tumor cells are c2 (subclone); remainder are c3.
# Shared SVs between c2+c3 = clonal; SVs private to each = subclonal.
c2_por=$(( ppur * sub_mix / 100 ))
c3_por=$(( ppur - c2_por ))
frac=($c2_por $c3_por $((100 - ppur)))

# --- VISOR seed ---
# Unique per (replicate, experiment, condition): no overlap
visor_seed=$(( (rep_idx + 1) * 10000 + (exp_idx + 1) * 1000 + (cond_idx + 1) * 10 ))

# --- machine-specific paths ---------------------------------------------------
# These paths were formerly embedded for one installation and did not resolve
# consistently on other machines. That concrete failure is why these scripts
# now consume the per-site configuration contract.
#
# They now come from config.local.sh, which is per-machine and never committed.
# The config is located by walking UP from this file, the one lookup that cannot
# itself be configured. $VISOR_CONFIG overrides it for SLURM tasks that have been
# copied to a node-local spool directory.
if [[ -z "${VISOR_CONFIG_LOADED:-}" ]]; then
    _cfg="${VISOR_CONFIG:-}"
    if [[ -z "$_cfg" ]]; then
        _d="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        while [[ "$_d" != "/" ]]; do
            [[ -r "$_d/config.local.sh" ]] && { _cfg="$_d/config.local.sh"; break; }
            _d="$(dirname "$_d")"
        done
    fi
    [[ -n "$_cfg" && -r "$_cfg" ]] || {
        echo "ERROR: config.local.sh not found by walking up from ${BASH_SOURCE[0]}." >&2
        echo "       cp config.example.sh config.local.sh at the repo root and edit it." >&2
        echo "       See README.md." >&2
        return 1 2>/dev/null || exit 1
    }
    # shellcheck disable=SC1090
    source "$_cfg"
    unset _cfg _d
fi

# --- Paths ---
ref="$REF_AUTO"
mc_hack_base="$MC_HACK_BASE"
hack_dir=${mc_hack_base}/${hack_name}
snp_hack=${mc_hack_base}/${snp_name}

# Normal BAM filenames differ between normal (sim.srt.bam) and o_normal (norm.bam)
norm_bam_files=(sim.srt.bam norm.bam norm.bam norm.bam norm.bam)
nbam=${VISOR_NBAM:-${NORM_SHORT_DIR}/${normal_name}/${norm_bam_files[$exp_idx]}}

# Replicate output. A downstream-only rerun can write into a fresh result root
# while reading large immutable inputs directly from a distinct source root.
# Keeping these roots separate removes the need for thousands of directory
# symlinks inside result trees.
base_dir="${VISOR_REPLICATES_DIR:-$REPLICATES_DIR}"
input_base_dir="${VISOR_INPUT_REPLICATES_DIR:-$base_dir}"
work_dir=${base_dir}/${rep}/${exp_name}/${cond_name}
input_work_dir=${input_base_dir}/${rep}/${exp_name}/${cond_name}
short_dir=${input_work_dir}/short
tbam=${VISOR_TBAM:-${short_dir}/short.out/sim.srt.bam}

manta_dir=${work_dir}/manta
svtyp_dir=${work_dir}/svtyp
snp_dir=${input_work_dir}/SNP
fac_dir=${input_work_dir}/facet
svc_dir=${work_dir}/svclone

# Shared scripts/tools. Workflow helpers must come from the versioned checkout,
# not from the output root: fresh REPLICATES_DIR roots intentionally contain
# data and results only.
samtools="${SAMTOOLS:-$(command -v samtools)}"
pipeline_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
helper_dir=${VISOR_REPLICATE_HELPER_DIR:-${pipeline_dir}/helpers}
get_sv_range=${helper_dir}/get_sv_range.R
make_input=${helper_dir}/make_input.R
svclone_cfg=${helper_dir}/svclone_config.ini
for helper in "$get_sv_range" "$make_input" "$svclone_cfg"; do
    [[ -r "$helper" ]] || {
        echo "ERROR: required VISOR helper is not readable: $helper" >&2
        return 1 2>/dev/null || exit 1
    }
done
seeded_py=$(dirname "$pipeline_dir")/tree_eval/eval_package/Phylogeny_benchmark/scripts/visor_seeded.py

# --- Reproducibility header ---
log_repro_header() {
    local script="${1:-unknown}"
    echo ""
    echo "========================================================"
    echo "  SCRIPT    : ${script}"
    echo "  DATE      : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  HOST      : $(hostname -s)"
    echo "  SLURM_JOB : ${SLURM_JOB_ID:-N/A}  TASK: ${SLURM_ARRAY_TASK_ID:-N/A}"
    echo "  rep       : ${rep}  (rep_idx=${rep_idx})"
    echo "  experiment: ${exp_name}  (exp_idx=${exp_idx})"
    echo "  condition : ${cond_name}  (cond_idx=${cond_idx})"
    echo "  purity    : ${ppur}%  mixture: ${sub_mix}%"
    echo "  clone frac: c2=${frac[0]}% c3=${frac[1]}% normal=${frac[2]}%"
    echo "  visor_seed: ${visor_seed}"
    echo "--------------------------------------------------------"
    echo "  hack_dir  : ${hack_dir}"
    echo "  snp_hack  : ${snp_hack}"
    echo "  nbam      : ${nbam}"
    echo "  tbam      : ${tbam}"
    echo "  input_dir : ${input_work_dir}"
    echo "  work_dir  : ${work_dir}"
    echo "========================================================"
    echo ""
}
