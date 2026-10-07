#!/bin/bash
# Submit VISOR replicate pipeline steps 01-03 with SLURM dependencies.
# 2250 tasks = 30 replicates x 5 experiments x 15 conditions
#             (p10/p20/p40/p60/p80 x m10/m30/m50)
#
# Step 00 (VISOR SHORtS) is NOT submitted here — all 2250 sim.srt.bam files
# are already generated (12 purity conditions via the original run; p10
# added separately via submit_10pct.sh). Do not re-run step 00.
#
# Dependency chain:
#   01 (Manta/SVtyper) → 02 (SNP pipeline)  [aftercorr = per-task]
#   02 (SNP pipeline)  → 03 (FACET)         [03 reads SNP.vcf.gz from 02]
# Step 04 (SVclone) is NOT submitted here.
# After 03 completes, run the SVclone arms and the SVCFit three-arm stages (see README.md).

set -euo pipefail

# --- machine-specific paths ---------------------------------------------------
# Absolute paths below come from config.local.sh, which is per-machine and never
# committed. Located by walking UP from this script. See README.md.
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
        exit 1
    }
    # shellcheck disable=SC1090
    source "$_cfg"
    unset _cfg _d
fi

cd "$(dirname "$0")"

mkdir -p log

# Verify all 2250 VISOR BAMs exist before submitting downstream steps.
n_bam=$(find "$REPLICATES_DIR"/rep* -name "sim.srt.bam" -path "*/short/short.out/sim.srt.bam" 2>/dev/null | wc -l)
if [ "$n_bam" -ne 2250 ]; then
    echo "ERROR: Found $n_bam/2250 sim.srt.bam files. Step 00 is not complete — cannot proceed."
    exit 1
fi
echo "Step 00 complete: $n_bam/2250 sim.srt.bam files present."
echo "Submitting steps 01-03 (2250 tasks each)..."
echo ""

JOB1=$(sbatch --parsable 01_manta_svtyper.sh)
echo "01_manta_svtyper   : $JOB1  (array 0-2249)"

# aftercorr: each task of JOB2 waits for its own corresponding task of JOB1,
# not for all 2250 Manta tasks. 02 needs svt_${cond_name}.vcf produced by 01.
JOB2=$(sbatch --parsable --dependency=aftercorr:${JOB1} 02_snp_pipeline.sh)
echo "02_snp_pipeline    : $JOB2  (array 0-2249, aftercorr $JOB1)"

JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} 03_facet.sh)
echo "03_facet           : $JOB3  (array 0-2249, afterok $JOB2 — needs SNP.vcf.gz)"

echo ""
echo "Monitor with: squeue -u \$USER"
echo ""
echo "After job 03 completes, run the SVclone arms and the SVCFit three-arm stages (see README.md)."
