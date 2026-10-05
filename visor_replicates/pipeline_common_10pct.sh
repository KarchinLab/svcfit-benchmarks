#!/bin/bash
set -euo pipefail
# Shared variables for p10 purity VISOR replicate pipeline.
# 30 replicates x 5 experiments x 3 conditions (p10m10, p10m30, p10m50)
# = 450 array tasks (0-449)
#
# Encoding:
#   rep_idx  = SLURM_ARRAY_TASK_ID / 15    (0-29 → rep1-rep30)
#   exp_idx  = (SLURM_ARRAY_TASK_ID % 15) / 3   (0-4 → exp1-exp5)
#   cond_idx = SLURM_ARRAY_TASK_ID % 3    (0-2 → p10m10, p10m30, p10m50)

# --- Decode task ID ---
rep_idx=$(( SLURM_ARRAY_TASK_ID / 15 ))
remainder=$(( SLURM_ARRAY_TASK_ID % 15 ))
exp_idx=$(( remainder / 3 ))
cond_idx=$(( remainder % 3 ))

rep=rep$(( rep_idx + 1 ))

# --- Experiment mapping (same as main pipeline) ---
exp_names=(exp1 exp2 exp3 exp4 exp5)
hack_names=(e1 e2 e3 e4 e5)
snp_names=(snp1 snp1_1 snp1_1 snp1_1 snp1_1)
normal_names=(normal o_normal o_normal o_normal o_normal)

exp_name=${exp_names[$exp_idx]}
hack_name=${hack_names[$exp_idx]}
snp_name=${snp_names[$exp_idx]}
normal_name=${normal_names[$exp_idx]}

# --- Condition mapping (p10 only) ---
mixtures=(10 30 50)

ppur=10
sub_mix=${mixtures[$cond_idx]}
cov=50
pur=100

# Condition name: c50p10m{sub_mix}
cond_name=c${cov}p${ppur}m${sub_mix}

# Clone fractions: same logic as main pipeline
c2_por=$(( ppur * sub_mix / 100 ))
c3_por=$(( ppur - c2_por ))
frac=($c2_por $c3_por $((100 - ppur)))

# --- VISOR seed ---
# Use cond offset +13..+15 to avoid collision with existing seeds (cond_idx+1 = 1-12 → *10 = 10-120)
visor_seed=$(( (rep_idx + 1) * 10000 + (exp_idx + 1) * 1000 + (cond_idx + 13) * 10 ))

# --- Paths (same as main pipeline) ---
# --- machine-specific paths ---------------------------------------------------
# The five paths below were hardcoded and did not resolve consistently. That is
# not hypothetical drift; it is the current state, and it is why these scripts
# could not be run here as committed.
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

ref="$REF_AUTO"
mc_hack_base="$MC_HACK_BASE"
hack_dir=${mc_hack_base}/${hack_name}
snp_hack=${mc_hack_base}/${snp_name}

norm_bam_files=(sim.srt.bam norm.bam norm.bam norm.bam norm.bam)
nbam=${NORM_SHORT_DIR}/${normal_name}/${norm_bam_files[$exp_idx]}

base_dir="$REPLICATES_DIR"
work_dir=${base_dir}/${rep}/${exp_name}/${cond_name}
short_dir=${work_dir}/short
tbam=${short_dir}/short.out/sim.srt.bam

manta_dir=${work_dir}/manta
svtyp_dir=${work_dir}/svtyp
snp_dir=${work_dir}/SNP
fac_dir=${work_dir}/facet
svc_dir=${work_dir}/svclone

samtools="${SAMTOOLS:-$(command -v samtools)}"
scrp_dir=${base_dir}/script
seeded_py=${base_dir}/visor_seeded.py

svclone_cfg=${base_dir}/script/svclone_config.ini
make_input=${base_dir}/script/make_input.R

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
    echo "  work_dir  : ${work_dir}"
    echo "========================================================"
    echo ""
}
