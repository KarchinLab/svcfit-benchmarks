#!/bin/bash
# Submit VISOR replicate pipeline steps 01-04 with SLURM dependencies.
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
#   01 + 03            → 04 (SVclone)
# After 03 completes, run the SVclone arms and the SVCFit three-arm stages (see README.md).

set -euo pipefail
cd "$(dirname "$0")"

mkdir -p log

# Verify all 2250 VISOR BAMs exist before submitting downstream steps.
n_bam=$(find rep* -name "sim.srt.bam" -path "*/short/short.out/sim.srt.bam" 2>/dev/null | wc -l)
if [ "$n_bam" -ne 2250 ]; then
    echo "ERROR: Found $n_bam/2250 sim.srt.bam files. Step 00 is not complete — cannot proceed."
    exit 1
fi
echo "Step 00 complete: $n_bam/2250 sim.srt.bam files present."
echo "Submitting steps 01-04 (2250 tasks each)..."
echo ""

JOB1=$(sbatch --parsable 01_manta_svtyper.sh)
echo "01_manta_svtyper   : $JOB1  (array 0-2249)"

# aftercorr: each task of JOB2 waits for its own corresponding task of JOB1,
# not for all 2250 Manta tasks. 02 needs svt_${cond_name}.vcf produced by 01.
JOB2=$(sbatch --parsable --dependency=aftercorr:${JOB1} 02_snp_pipeline.sh)
echo "02_snp_pipeline    : $JOB2  (array 0-2249, aftercorr $JOB1)"

JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} 03_facet.sh)
echo "03_facet           : $JOB3  (array 0-2249, afterok $JOB2 — needs SNP.vcf.gz)"

echo "04_svclone         : $JOB4  (array 0-2249, afterok $JOB1+$JOB3)"

echo ""
echo "Monitor with: squeue -u \$USER"
echo ""
echo "After jobs 02 + 04 complete, run the R analysis:"
echo "  then run the SVclone arms and the SVCFit three-arm stages described in README.md"
