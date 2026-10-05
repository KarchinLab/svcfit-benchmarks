#!/bin/bash
#SBATCH --job-name=longi_short
#SBATCH --output=log/longi_short_%A_%a.out
#SBATCH --error=log/longi_short_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --array=0-4

###############################################################################
# longi_short.sh
#
# Generates TWO simulated BAMs per scenario (timepoint 1 and timepoint 2)
# using VISOR SHORtS, with different clone fractions at each timepoint.
#
# This reuses the existing HACk genomes from the single-timepoint benchmark.
# The only difference is that each scenario now produces two BAMs with
# different --clonefraction values to simulate longitudinal sampling.
#
# Usage:
#   Set PPUR, SCENARIO, and COV below, then submit:
#     sbatch longi_short.sh
#
# The SLURM array iterates over exp1..exp5 (SV-CNV overlap configurations).
#
# S2 is a negative control (identical fractions at both timepoints); in S1, S3 and S4 the
# two subclones have distinct CCF trajectories.
###############################################################################

thread=$SLURM_CPUS_PER_TASK
cov=50

########## DEFAULTS (override via --export from run_all.sh) ##########
# Tumor purity (percent of cells that are tumor)
ppur=${ppur:-60}

# Longitudinal scenario index (S1-S4), controls clone dynamics
# S1: large subclone swap         — c2 dominant (5:1) swaps to c3 dominant (1:5)
# S2: negative control            — identical fractions at both timepoints (no clonal evolution)
# S3: gradual drift               — equal subclones at t1; c2 expands at t2 (purity constant)
# S4: asymmetric regression       — equal subclones at t1; purity halves, c2 resists better (2:1 at t2)
SCENARIO=${SCENARIO:-S1}

# Bootstrap run index.  BOOT=0 uses the original seeds and directory structure
# (backward compatible with existing data).  BOOT=1..5 produce new replicates
# in b{BOOT}/ subdirectories with distinct seeds.
BOOT=${BOOT:-0}
######################################################################

# Clone fractions for c2, c3, and normal (must sum to 100)
# These represent the fraction of cells in the tumor BAM.
# c2 and c3 each contain truncal (c1) SVs plus their private SVs.
# The "normal" fraction has only germline SNPs (no tumor SVs).
case $SCENARIO in
  S1)
    # Subclone swap: c2 dominant (5:1) -> c3 dominant (1:5)
    T1_C2=$((ppur * 5 / 6)); T1_C3=$((ppur - T1_C2))
    T2_C2=$((ppur * 1 / 6)); T2_C3=$((ppur - T2_C2))
    ;;
  S2)
    # Negative control: identical fractions at both timepoints — no clonal evolution.
    # DP-GMM space: sub1 and sub2 occupy the same point (ppur/2, ppur/2); expected to merge.
    T1_C2=$((ppur / 2)); T1_C3=$((ppur - T1_C2))
    T2_C2=$((ppur / 2)); T2_C3=$((ppur - T2_C2))
    ;;
  S3)
    # Gradual clonal drift: equal subclones at t1; c2 expands strongly, c3 contracts.
    # Total purity unchanged between timepoints.
    # t1: c2=1/2*ppur, c3=1/2*ppur   t2: c2=4/5*ppur, c3=1/5*ppur
    # DP-GMM space: sub1=(+3/10, ppur/2), sub2=(-3/10, ppur/2) — same t1, delta=±0.30
    T1_C2=$((ppur / 2));     T1_C3=$((ppur - T1_C2))
    T2_C2=$((ppur * 4 / 5)); T2_C3=$((ppur - T2_C2))
    ;;
  S4)
    # Asymmetric regression: equal subclones at t1; purity halves at t2 but c2 resists
    # better (c2:c3 = 4:1 within tumor at t2 vs 1:1 at t1).
    # t1: c2=1/2*ppur, c3=1/2*ppur   t2: c2=2/5*ppur, c3=1/10*ppur
    # DP-GMM space: sub1=(+3/10, ppur/2), sub2=(-3/10, ppur/2) — delta=±0.30, same as S3
    # T2_C2+T2_C3 = ppur/2 (purity halves); CCF within tumor = 4/5 vs 1/5.
    T1_C2=$((ppur / 2));     T1_C3=$((ppur - T1_C2))
    T2_C2=$((ppur * 2 / 5)); T2_C3=$((ppur / 10))
    ;;
  *)
    echo "Unknown scenario: $SCENARIO" && exit 1
    ;;
esac

# Compute normal fractions from tumor clone fractions (must sum to 100)
T1_NORM=$((100 - T1_C2 - T1_C3))
T2_NORM=$((100 - T2_C2 - T2_C3))
frac_t1=($T1_C2 $T1_C3 $T1_NORM)
frac_t2=($T2_C2 $T2_C3 $T2_NORM)

# Experiment (SV-CNV overlap config)
exp_lst=(exp1 exp2 exp3 exp4 exp5)
hack_lst=(e1 e2 e3 e4 e5)
snp_lst=(snp1 snp1_1 snp1_1 snp1_1 snp1_1)
nor_lst=(normal o_normal o_normal o_normal o_normal)
normal_name=${nor_lst[$SLURM_ARRAY_TASK_ID]}
samp_name=${exp_lst[$SLURM_ARRAY_TASK_ID]}
hack_file=${hack_lst[$SLURM_ARRAY_TASK_ID]}

# --- PATHS ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PUB_DIR="$(dirname "$SCRIPT_DIR")"
hack_base=$PUB_DIR/data                                       # hack/ lives under publishable/data/
out_base=$PUB_DIR/outputs                                     # simulation output root
ref=$PUB_DIR/data/reference/chr1-2.fa
seeded_py=$SCRIPT_DIR/visor_seeded.py

snp_hack=$hack_base/hack/${snp_lst[$SLURM_ARRAY_TASK_ID]}
hack_dir=$hack_base/hack/$hack_file

# BOOT=0: original directory (backward compatible); BOOT>0: b{BOOT}/ subdirectory
if [[ $BOOT -eq 0 ]]; then
  longi_dir=$out_base/${SCENARIO}/c${cov}p${ppur}/$samp_name
else
  longi_dir=$out_base/${SCENARIO}/c${cov}p${ppur}/b${BOOT}/$samp_name
fi

# --- Deterministic seeds (design v4) ---
# Bootstrap extension: BOOT=0 uses the original formula (seeds 500000-543052)
#   so that existing BAMs are bit-for-bit reproducible.
# BOOT=1..5 add a 100000*BOOT offset, placing each bootstrap run in a
#   non-overlapping seed range (600000-643052, 700000-743052, …, 1000000-1043052).
# No collision with visor_replicates (≤315150), v1 seeds (400000-449952),
#   or any prior range.
# formula: base = 500000 + BOOT*100000 + scen_idx*10000 + pur_tier*1000 + (exp_idx+1)*10
#          T1 = base+1, T2 = base+2
scen_idx=$(echo "$SCENARIO" | tr -d 'S')    # S1→1 ... S4→4
case $ppur in 40) pur_tier=1 ;; 60) pur_tier=2 ;; 80) pur_tier=3 ;; *) pur_tier=0 ;; esac
exp_idx=$SLURM_ARRAY_TASK_ID
base_seed=$(( 500000 + BOOT * 100000 + scen_idx * 10000 + pur_tier * 1000 + (exp_idx + 1) * 10 ))
seed_t1=$(( base_seed + 1 ))
seed_t2=$(( base_seed + 2 ))

REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
# shellcheck disable=SC1091
source "${VISOR_CONFIG:-$REPO_ROOT/config.local.sh}"
source "$CONDA_SH"
conda activate "$ENV_VISOR"

mkdir -p $longi_dir/t1 $longi_dir/t2
mkdir -p log

###############################################################################
# Timepoint 1
###############################################################################
echo "=== Timepoint 1: clonefraction = ${frac_t1[@]} ==="

cut -f1,2 \
    $hack_dir/c2/*.fai $hack_dir/c3/*.fai $snp_hack/*.fai \
    | sort \
    | awk '$2 > maxvals[$1] {lines[$1]=$0; maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' \
    | awk 'OFS=FS="\t"{print $1, "1", $2, 100, 100}' > $longi_dir/t1/short.bed

rm -rf $longi_dir/t1/short.out
echo "Running VISOR SHORtS T1 with seed=${seed_t1}..."
"$ENV_VISOR/bin/python" $seeded_py $seed_t1 \
    -g $ref \
    -s $hack_dir/c2 $hack_dir/c3 $snp_hack \
    -b $longi_dir/t1/short.bed \
    -o $longi_dir/t1/short.out \
    --mutation 0 \
    --extindels 0 \
    --threads $thread \
    --coverage $cov \
    --clonefraction ${frac_t1[@]} \
    --error 0 \
    --indels 0

###############################################################################
# Timepoint 2
###############################################################################
echo "=== Timepoint 2: clonefraction = ${frac_t2[@]} ==="

cut -f1,2 \
    $hack_dir/c2/*.fai $hack_dir/c3/*.fai $snp_hack/*.fai \
    | sort \
    | awk '$2 > maxvals[$1] {lines[$1]=$0; maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' \
    | awk 'OFS=FS="\t"{print $1, "1", $2, 100, 100}' > $longi_dir/t2/short.bed

rm -rf $longi_dir/t2/short.out
echo "Running VISOR SHORtS T2 with seed=${seed_t2}..."
"$ENV_VISOR/bin/python" $seeded_py $seed_t2 \
    -g $ref \
    -s $hack_dir/c2 $hack_dir/c3 $snp_hack \
    -b $longi_dir/t2/short.bed \
    -o $longi_dir/t2/short.out \
    --mutation 0 \
    --extindels 0 \
    --threads $thread \
    --coverage $cov \
    --clonefraction ${frac_t2[@]} \
    --error 0 \
    --indels 0

t1bam=$longi_dir/t1/short.out/sim.srt.bam
t2bam=$longi_dir/t2/short.out/sim.srt.bam
[[ ! -f $t1bam ]] && { echo "ERROR: T1 BAM not produced" >&2; exit 1; }
[[ ! -f $t2bam ]] && { echo "ERROR: T2 BAM not produced" >&2; exit 1; }

echo "=== Done: $SCENARIO / $samp_name ==="
echo "T1 BAM: $t1bam  ($(du -sh $t1bam 2>/dev/null | cut -f1))"
echo "T2 BAM: $t2bam  ($(du -sh $t2bam 2>/dev/null | cut -f1))"
