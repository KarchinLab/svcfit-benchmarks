#!/bin/bash
# Submit N replicates of the c50 simulation, in waves.
#
#   COV=50 ./16_chrx_replicates_all.sh [N] [WAVE] [FIRST]
#     N     how many replicates          (default 30)
#     WAVE  replicates started at once   (default 6)
#     FIRST first replicate index        (default 1)
#
# WHY WAVES. The read simulation is 45 tasks x 12 cpu = 540 cpu per replicate,
# against a `normal` QoS ceiling of cpu=3600. Six replicates is 3,240 -- under
# the cap with room for the 1-cpu segmentation tasks to interleave. Submitting
# all 30 at once would not break anything (SLURM would just hold them PENDING on
# QOSMaxCpuPerUserLimit) but it puts ~4,100 tasks in a shared queue at once and
# makes the run impossible to read in squeue.
#
# Each wave's read simulation waits on the previous wave's with `afterany`, so a
# failed replicate stalls nothing behind it. Downstream stages of earlier waves
# overlap later waves' simulation, which is where most of the wall clock is
# recovered.
#
# IDEMPOTENT, per stage rather than per replicate: 02, 12 and 06 each skip a task
# whose output already exists, so a wave that half-failed can be resubmitted.

set -euo pipefail

N="${1:-30}"
WAVE="${2:-6}"
FIRST="${3:-1}"
COV="${COV:?set COV, e.g. COV=50 $0}"

for v in N WAVE FIRST; do
    [[ "${!v}" =~ ^[0-9]+$ ]] || { echo "ERROR: $v must be an integer, got '${!v}'" >&2; exit 1; }
done
(( FIRST > 0 )) || { echo "ERROR: replicates start at 1 (REP=0 is the original run)" >&2; exit 1; }
(( WAVE > 0 )) || { echo "ERROR: WAVE must be positive" >&2; exit 1; }

_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_DIR}/chrx_common.sh"

last=$(( FIRST + N - 1 ))
echo "replicates ${FIRST}..${last} at COV=${COV}, ${WAVE} per wave"
echo

prev=""          # the previous wave's read-simulation job, for the throttle
i=0
for (( rep = FIRST; rep <= last; rep++ )); do
    if (( i % WAVE == 0 )) && [[ -n "$prev" ]]; then
        echo "--- wave boundary: holding on ${prev} ---"
        hold="$prev"
    elif (( i % WAVE == 0 )); then
        hold=""
    fi
    COV="$COV" "${_DIR}/15_chrx_replicate_chain.sh" "$rep" ${hold:+"$hold"}
    # The last chain submitted in this wave becomes the next wave's gate.
    (( (i + 1) % WAVE == 0 )) && prev="$(cat "${BASE_DIR}/log/.rep${rep}.shorts_jobid")"
    i=$(( i + 1 ))
done

echo
echo "submitted ${N} replicate chains. Watch with:"
echo "  squeue -u \$(whoami) -o '%.14i %.12j %.9T %.10M %R'"
