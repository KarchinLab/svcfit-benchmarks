#!/bin/bash
#
# Submit the extract array, then submit the aggregator with afterok dependency.
# Usage: bash run_all.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"

EXTRACT_JOB=$(sbatch --parsable "$SCRIPT_DIR/submit_extract_array.sh")
echo "Submitted extract array : $EXTRACT_JOB"

AGG_JOB=$(sbatch --parsable --dependency=afterok:$EXTRACT_JOB "$SCRIPT_DIR/submit_aggregate.sh")
echo "Submitted aggregate     : $AGG_JOB  (depends on $EXTRACT_JOB)"

cat <<EOF

Monitor with:   squeue -u \$USER
Logs:           $SCRIPT_DIR/logs/
Outputs:        will appear under \$OUT_DIR (see config.R)

When the aggregate job finishes, the headline numbers are in:
  - a2_spearman_by_clone.csv      (3 rho values + 95% CIs)
  - a2_lmer_summary.txt           (mixed model + interaction LRT p-value)
  - a2_interpretation.txt         (one-line decision-rule verdict)
  - fig_a2_3.pdf                  (per-clone scatter)

EOF
