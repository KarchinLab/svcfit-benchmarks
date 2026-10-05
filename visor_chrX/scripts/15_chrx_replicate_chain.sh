#!/bin/bash
# Submit the whole pipeline for ONE replicate, as a dependency chain.
#
#   ./15_chrx_replicate_chain.sh <REP> [AFTER_JOBID]
#
# REP selects an independent draw of SHORtS read noise -- the quantity a block
# bootstrap has to resample. REP=0 is the original run and is refused here: it
# already exists, and re-simulating it would overwrite the committed c50 results.
#
# AFTER_JOBID, if given, holds this replicate's read simulation until that job
# finishes (afterany, not afterok -- a wave is a throttle, not a prerequisite).
# 16_chrx_replicates_all.sh uses it to keep concurrent SHORtS within the QoS cap.
#
# Stages, and why the dependencies are shaped this way:
#
#   A  02  read simulation, 45 tasks x 12 cpu        <- the seeded stage
#   B  12  Manta + SVtyper, 45 tasks   after A       \ independent of each other,
#   C  06  depth segmentation, 45 tasks after A      / so both hang off A
#   D  17  join + score, 1 task        after B and C
#
# B and C are deliberately siblings rather than a sequence: nothing in the
# segmentation reads the calls, and serialising them would add an hour per
# replicate for nothing.
#
# THE MATCHED NORMAL IS SHARED, NOT RESIMULATED. 02's normal mode uses a
# hardcoded seed (999) that REP does not reach, so a per-replicate normal would
# be byte-identical to replicate 0's anyway unless that seed is also made
# REP-dependent -- which was considered and deliberately not done. The
# consequence belongs in any write-up: cn_bar is a tumour-to-normal depth ratio,
# so with one frozen normal the bootstrap resamples TUMOUR read noise only, and
# intervals on the cn_bar-dependent branches (e2, e4, depth-routed deletions)
# will be narrower than a full resampling would give.

set -euo pipefail

REP="${1:?usage: $0 <REP> [AFTER_JOBID]}"
AFTER="${2:-}"

[[ "$REP" =~ ^[0-9]+$ ]] || { echo "ERROR: REP must be a non-negative integer, got '$REP'" >&2; exit 1; }
(( REP > 0 )) || { echo "ERROR: REP=0 is the original run and already exists; refusing to overwrite it." >&2
                   echo "       Replicates start at 1." >&2; exit 1; }

COV="${COV:?set COV, e.g. COV=50}"

_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load the machine configuration before choosing the run root. chrx_common.sh
# defaults BASE_DIR to the configured CHRX_DIR.
source "${_DIR}/chrx_common.sh"
cd "$BASE_DIR"                      # log/ in the SBATCH directives is relative
mkdir -p log

# For CHRX_EXCLUDE_NODES -- the nodes that kill tasks at start-up. Defined in
# chrx_common.sh rather than here so 15 and 18 cannot drift apart on it.
excl=""; [[ -n "${CHRX_EXCLUDE_NODES:-}" ]] && excl="--exclude=${CHRX_EXCLUDE_NODES}"

short_rep="${BASE_DIR}/short_rep${REP}"
shared_normal="${BASE_DIR}/short/normal_c${COV}"

# The shared normal, reached at the path each replicate's scripts expect. 06 and
# 11 both look for ${short_base}/normal_c${cov}, and short_base carries the
# replicate tag, so without this the replicate has no normal at all.
[[ -d "$shared_normal" ]] || { echo "ERROR: no shared normal at $shared_normal" >&2; exit 1; }
mkdir -p "$short_rep"
if [[ ! -e "${short_rep}/normal_c${COV}" ]]; then
    ln -s "$shared_normal" "${short_rep}/normal_c${COV}"
    echo "  linked normal_c${COV} -> ${shared_normal}"
fi

EXP="ALL,COV=${COV},REP=${REP}"
dep_a=""; [[ -n "$AFTER" ]] && dep_a="--dependency=afterany:${AFTER}"

# A: read simulation. The only stage the seed reaches.
a=$(sbatch --parsable $dep_a $excl --export="$EXP" --array=0-44%45 \
        scripts/02_chrx_shorts.sh)

# B: calling.  C: segmentation.  Both afterok:A -- there is nothing to call or
# segment if the reads did not simulate.
b=$(sbatch --parsable --dependency=afterok:$a $excl --export="$EXP" --array=0-44%45 \
        scripts/12_chrx_calls_array.sbatch)
c=$(sbatch --parsable --dependency=afterok:$a $excl \
        --export="${EXP},BIN_KB=10,SEG_DIR=${BASE_DIR}/seg10k" \
        --time=6:00:00 --array=0-44%45 \
        scripts/06_chrx_segment_array.sbatch)

# D: join + score, once both are in.
d=$(sbatch --parsable --dependency=afterok:${b}:${c} $excl --export="$EXP" \
        scripts/17_chrx_score_one.sbatch)

printf 'rep %-3s  shorts=%s  calls=%s  seg=%s  score=%s\n' "$REP" "$a" "$b" "$c" "$d"
printf '%s\n' "$a" > "${BASE_DIR}/log/.rep${REP}.shorts_jobid"
