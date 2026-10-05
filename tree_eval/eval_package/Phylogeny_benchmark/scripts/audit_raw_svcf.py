#!/usr/bin/env python3
"""
Audit of `mean raw_svcf / p` for non-DUP trunk SVs in diploid background segments.

Purpose
-------
Verify (and, if necessary, correct) the empirical values in
`longitudinal/evaluation/boots0_4/ccf_bias_assessment.md` §"Step 1"
(and the identical table in `publishability_assessment.md` §"Part 3, Step 1").

The original audit reported:
    p=20%: n=277, mean/p = 1.014
    p=40%: n=252, mean/p = 0.880
    p=60%: n=263, mean/p = 0.899
    p=80%: n=268, mean/p = 0.909
but the script that produced those numbers is not present in the tree; only
the summary table is stored. This script reproduces the computation from
first principles using the actual BED files and the VISOR simulation
ground truth for trunk SV positions.

Method
------
For each purity p ∈ {20, 40, 60, 80}%:
  1. Identify ground-truth trunk SV positions from the VISOR input files
     `VISOR_longitudinal/publishable/data/hack/{c1,c11}.bed` (both trunk
     haplotypes; 66 unique positions total).
  2. For each of the 25 replicates (5 bootstraps × 5 experiments) at that
     purity, iterate over both timepoint BED files
     (`svcfit_output/S1_p{P}_exp{E}_t{T}.bed`).
  3. Filter each BED file to rows where:
       - `classification != 'DUP'`     (excluding duplications; formula differs)
       - `bkg_cnv == 'norm'`           (diploid background segment)
       - `CHROM, POS` match a truth-trunk position within POS_TOL bp
  4. Deduplicate SVs within each BED file by `(CHROM, POS, END)` because
     Manta records BND / paired breakpoints in two rows sharing the same
     raw_svcf value.
  5. Compute `mean(raw_svcf)` pooled across all matched SVs at that purity,
     and divide by p_frac = p/100 to obtain the ratio.

The `raw_svcf` column in the BED is the value produced by SVCFit's
`calc_svcf.R`:
    raw_svcf = 2 × sv_alt / (sv_alt + sv_ref)    if classification != 'DUP'
    raw_svcf = raw_svcf / 2                       if zygosity == 'hom'
No further correction is applied here (i.e., the halving of hom-called SVs
is retained). For a truly heterozygous trunk SV at purity p, the expected
value is p; for a truly homozygous trunk SV called as hom, also p; for a
truly het SV mis-called as hom by SVtyper, ~p/2. Systematic zygosity
misclassification therefore biases the ratio downward.

Outputs
-------
CSV: audit_raw_svcf_results.csv
    Columns: purity, n_svs, mean_raw_svcf, ratio_raw_svcf_over_p,
             n_replicates, n_bed_files
"""

import csv
import os
import statistics
import sys
from collections import defaultdict

# Paths resolve relative to this package so it runs anywhere. Override with
# SVCFIT_EVAL_BASE / SVCFIT_EVAL_TRUTH / SVCFIT_EVAL_OUT to point elsewhere.
_PKG       = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_DIR   = os.environ.get('SVCFIT_EVAL_BASE',
                            os.path.join(_PKG, 'bootstrap', 'S1'))
TRUTH_DIR  = os.environ.get('SVCFIT_EVAL_TRUTH',
                            os.path.join(_PKG, 'input_data', 'hack'))
OUT_PATH   = os.environ.get('SVCFIT_EVAL_OUT',
                            os.path.join(_PKG, 'output', 'boots0_4',
                                         'audit_raw_svcf_results.csv'))
PURITIES   = [20, 40, 60, 80]
# "boots0_4" means boots 0..4 = 5 bootstraps.  Boot 0 is stored at the top level
# (c50p{P}/exp{E}/), boots 1..4 at c50p{P}/b{B}/exp{E}/.
BOOTS      = [0, 1, 2, 3, 4]
EXPS       = list(range(1, 6))          # exp1..exp5
TIMEPOINTS = [1, 2]                     # t1, t2
POS_TOL    = 5                          # bp tolerance for truth-vs-called POS


def bed_path(purity, boot, exp, tp):
    """Boot 0 lives at c50p{P}/exp{E}/; boots >=1 live at c50p{P}/b{B}/exp{E}/."""
    if boot == 0:
        return (f'{BASE_DIR}/c50p{purity}/exp{exp}/svcfit_output/'
                f'S1_p{purity}_exp{exp}_t{tp}.bed')
    return (f'{BASE_DIR}/c50p{purity}/b{boot}/exp{exp}/svcfit_output/'
            f'S1_p{purity}_exp{exp}_t{tp}.bed')


def load_trunk_positions():
    """Ground-truth trunk SV start positions from VISOR hack files."""
    pos = set()
    for hap in ('c1', 'c11'):
        with open(os.path.join(TRUTH_DIR, f'{hap}.bed')) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 4:
                    pos.add((parts[0], int(parts[1])))
    return pos


def is_trunk(chrom, pos, trunk_set, tol=POS_TOL):
    for c, p in trunk_set:
        if c == chrom and abs(p - pos) <= tol:
            return True
    return False


def audit_one_purity(purity, trunk_set):
    p_frac = purity / 100.0
    values = []
    n_files = 0
    reps_seen = set()
    for b in BOOTS:
        for e in EXPS:
            for tp in TIMEPOINTS:
                bed = bed_path(purity, b, e, tp)
                if not os.path.exists(bed):
                    continue
                n_files += 1
                reps_seen.add((b, e))
                seen = {}
                with open(bed) as fh:
                    for row in csv.DictReader(fh, delimiter='\t'):
                        if row.get('classification') == 'DUP':
                            continue
                        if row.get('bkg_cnv') != 'norm':
                            continue
                        try:
                            pos = int(row['POS'])
                            raw = float(row['raw_svcf'])
                        except (ValueError, KeyError, TypeError):
                            continue
                        if not is_trunk(row['CHROM'], pos, trunk_set):
                            continue
                        key = (row['CHROM'], pos, row.get('END', '?'))
                        if key not in seen:
                            seen[key] = raw
                values.extend(seen.values())
    if not values:
        return None
    mean_raw = sum(values) / len(values)
    median_raw = statistics.median(values)
    # 10% trimmed mean
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    lo, hi = int(n * 0.10), n - int(n * 0.10)
    trimmed = sorted_vals[lo:hi]
    trim_mean = sum(trimmed) / len(trimmed) if trimmed else float('nan')
    return {
        'purity': purity,
        'n_svs': len(values),
        'mean_raw_svcf': round(mean_raw, 4),
        'ratio_mean_over_p': round(mean_raw / p_frac, 4),
        'ratio_median_over_p': round(median_raw / p_frac, 4),
        'ratio_trimmean10_over_p': round(trim_mean / p_frac, 4),
        'n_replicates': len(reps_seen),
        'n_bed_files': n_files,
    }


def main():
    trunk_set = load_trunk_positions()
    print(f'Loaded {len(trunk_set)} ground-truth trunk positions '
          f'from {TRUTH_DIR}/{{c1,c11}}.bed', file=sys.stderr)

    results = []
    for p in PURITIES:
        r = audit_one_purity(p, trunk_set)
        if r is not None:
            results.append(r)
            print(f'p={p:2d}%: n={r["n_svs"]:4d} SVs across '
                  f'{r["n_replicates"]} replicates '
                  f'({r["n_bed_files"]} BED files); '
                  f'ratio_mean/p={r["ratio_mean_over_p"]:.3f}, '
                  f'median/p={r["ratio_median_over_p"]:.3f}, '
                  f'trim10%/p={r["ratio_trimmean10_over_p"]:.3f}',
                  file=sys.stderr)

    with open(OUT_PATH, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f'\nWrote {OUT_PATH}', file=sys.stderr)


if __name__ == '__main__':
    main()
