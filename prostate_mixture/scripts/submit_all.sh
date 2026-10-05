#!/bin/bash
# Submit full 30-replicate prostate mixture pipeline with SLURM dependencies.
# Dependency chain:
#   00a (filter source BAMs, run once) → 00 (make mixtures)
#   00  → 01 (manta/svtyper), 02 (GATK SNP call)
#   01 + 02 → 03 (SNP process)
#   03 → 04 (FACET; needs merged SNP VCF from 03)
#   01 + 04 → 05 (SVclone)
# After 05 completes, run 04a/04b (chrX), submit_svcfit_correction.sh and 08_prostate_three_arm.R (see README.md).
# Each step script is idempotent: tasks whose output already exists exit immediately,
# so resubmitting the full pipeline is safe and will skip completed work.

set -euo pipefail
cd "$(dirname "$0")"

mkdir -p log

echo "Submitting 30-replicate prostate mixture pipeline..."

# pipeline_common.sh is sourced here only for src_bams, so the array decoding is not meaningful yet.
: "${SLURM_ARRAY_TASK_ID:=0}"
# shellcheck disable=SC1091
source ./pipeline_common.sh
if [[ -f ${src_bams}/fil_bM.bam && -f ${src_bams}/fil_gM.bam && \
      -f ${src_bams}/fil_bM_odd.bam && -f ${src_bams}/fil_bM_even.bam && \
      -f ${src_bams}/fil_gM_odd.bam && -f ${src_bams}/fil_gM_odd_downs.bam && \
      -f ${src_bams}/fil_gM_even_downs.bam ]]; then
    echo "00a: all filtered source BAMs already exist, skipping."
    DEP_00A=""
else
    JOB0A=$(sbatch --parsable 00a_filter_source_bams.sh)
    echo "00a_filter_source  : $JOB0A  (array 0-6, ~16h for 130G BAMs)"
    DEP_00A="--dependency=afterok:${JOB0A}"
fi

JOB0=$(sbatch --parsable $DEP_00A 00_make_mixtures.sh)
echo "00_make_mixtures   : $JOB0  (array 0-329, after ${JOB0A:-skipped})"

JOB1=$(sbatch --parsable --dependency=afterok:${JOB0} 01_manta_svtyper.sh)
echo "01_manta_svtyper   : $JOB1  (array 0-329, after $JOB0)"

JOB2=$(sbatch --parsable --dependency=afterok:${JOB0} 02_snp_call.sh)
echo "02_snp_call        : $JOB2  (array 0-6929, after $JOB0)"

JOB3=$(sbatch --parsable --dependency=afterok:${JOB1},afterok:${JOB2} 03_snp_process.sh)
echo "03_snp_process     : $JOB3  (array 0-329, after $JOB1 + $JOB2)"

JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} 04_facet.sh)
echo "04_facet           : $JOB4  (array 0-329, after $JOB3)"

JOB5=$(sbatch --parsable --dependency=afterok:${JOB1},afterok:${JOB4} 05_svclone.sh)
echo "05_svclone         : $JOB5  (array 0-329, after $JOB1 + $JOB4)"

echo ""
echo "All jobs submitted. Monitor with:"
echo "  squeue -u \$USER"
echo ""
echo "After job $JOB5 completes, run the chrX stages (04a, 04b), submit_svcfit_correction.sh"
echo "and 08_prostate_three_arm.R, as described in README.md."
