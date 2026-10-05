#!/usr/bin/env python3
"""Score the chrX segmentation against analytic truth, across the condition grid.

The simulation is the only place the depth segmentation can be checked
against a known answer: on the real cohort the true copy number is what we are
trying to measure. Here it is fixed by construction.

Normal chrX is single-copy, so for a locus L:

    cn_bar(L) = f_c2 * cn_c2(L) + f_c3 * cn_c3(L) + f_normal * 1

with cn_c2 / cn_c3 read straight out of the amplification beds
(`truth/sv_beds/dup_c22.bed`, `dup_c33.bed`) and 1 outside them. `info = N` in a
VISOR `tandem duplication` row means N TOTAL copies -- HACk.py:197 emits the
segment `inf` times, replacing the original -- so on hemizygous chrX the local
copy number is N directly.

e1 carries no amplification at all, so its predicted cn_bar is 1.00 everywhere.

REGION INTERIORS ONLY. Segment boundaries land on 100 kb bin edges, so the bins
straddling an amplification edge are mixtures. Each region is trimmed by
--trim (default 500 kb) at both ends before the observed level is measured;
otherwise the edge bins drag the mean and it looks like a bias when it is just
resolution.

e2 IS EXPECTED TO READ SLIGHTLY BELOW e4 at matched conditions, and below
prediction in both. In e2 the SVs precede the amplification, so a deletion is
absent from every copy (cn 0 there); in e4 the amplification comes first and a
deletion removes only one copy of a copies. Neither is modelled here -- the
prediction is the amplification level alone. The prior session measured this
same effect in the raw depth ratios: e2 below e4 in 15 of 15 conditions.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

EXPERIMENT_HAS_AMP = {"e1": False, "e2": True, "e4": True}


def read_dup_bed(path):
    """(start, end) -> total copies, from a VISOR tandem-duplication bed."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[3] != "tandem duplication":
                continue
            out[(int(f[1]), int(f[2]))] = int(f[4])
    return out


def read_seg(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            rows.append((int(r["start"]), int(r["end"]), float(r["cn_bar"])))
    return rows


def observed_level(seg, lo, hi):
    """Length-weighted mean cn_bar over segments overlapping [lo, hi)."""
    num = den = 0.0
    for s, e, cn in seg:
        ov = min(e, hi) - max(s, lo)
        if ov > 0:
            num += cn * ov
            den += ov
    return (num / den) if den > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditions", required=True)
    ap.add_argument("--sv-bed-dir", required=True)
    ap.add_argument("--seg-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--trim", type=int, default=500_000,
                    help="bp trimmed from each region end before measuring")
    args = ap.parse_args()

    c2_amp = read_dup_bed(os.path.join(args.sv_bed_dir, "dup_c22.bed"))
    c3_amp = read_dup_bed(os.path.join(args.sv_bed_dir, "dup_c33.bed"))

    # Union of amplified intervals; copies default to 1 in a clone that lacks it.
    regions = []
    for iv in sorted(set(c2_amp) | set(c3_amp)):
        regions.append({"start": iv[0], "end": iv[1],
                        "cn_c2": c2_amp.get(iv, 1), "cn_c3": c3_amp.get(iv, 1)})
    label = {}
    for r in regions:
        if r["cn_c2"] > 1 and r["cn_c3"] > 1:
            label[(r["start"], r["end"])] = "A_clonal"
        elif r["cn_c2"] > 1:
            label[(r["start"], r["end"])] = "B_c2_private"
        else:
            label[(r["start"], r["end"])] = "C_c3_private"
    print(f"amplification regions from the dup beds ({len(regions)}):")
    for r in regions:
        print(f"  {label[(r['start'], r['end'])]:<14} "
              f"{r['start']:>11,}-{r['end']:>11,}  cn_c2={r['cn_c2']}  cn_c3={r['cn_c3']}")

    with open(args.conditions) as fh:
        conditions = list(csv.DictReader(fh, delimiter="\t"))

    rows = []
    missing = []
    for cond in sorted(conditions, key=lambda c: int(c["task_id"])):
        exp, cname = cond["experiment"], cond["condition"]
        seg_path = os.path.join(args.seg_dir, f"chrx_seg_{exp}_{cname}.csv")
        if not os.path.exists(seg_path):
            missing.append(f"{exp}/{cname}")
            continue
        seg = read_seg(seg_path)
        f2 = int(cond["c2_pct"]) / 100
        f3 = int(cond["c3_pct"]) / 100
        fn = int(cond["normal_pct"]) / 100

        targets = [("baseline", None, 1, 1)]
        if EXPERIMENT_HAS_AMP[exp]:
            for r in regions:
                targets.append((label[(r["start"], r["end"])],
                                (r["start"], r["end"]), r["cn_c2"], r["cn_c3"]))

        for name, iv, cn2, cn3 in targets:
            pred = f2 * cn2 + f3 * cn3 + fn * 1
            if iv is None:
                # baseline = the stretches outside every amplified region
                if EXPERIMENT_HAS_AMP[exp]:
                    gaps, prev = [], 0
                    for r in regions:
                        if r["start"] > prev:
                            gaps.append((prev, r["start"]))
                        prev = max(prev, r["end"])
                    gaps.append((prev, 156_040_895))
                    gaps = [(a + args.trim, b - args.trim) for a, b in gaps
                            if b - a > 2 * args.trim]
                else:
                    gaps = [(args.trim, 156_040_895 - args.trim)]
                num = den = 0.0
                for a, b in gaps:
                    v = observed_level(seg, a, b)
                    if v == v:
                        num += v * (b - a)
                        den += (b - a)
                obs = num / den if den else float("nan")
            else:
                obs = observed_level(seg, iv[0] + args.trim, iv[1] - args.trim)

            rows.append({
                "experiment": exp, "condition": cname, "task_id": cond["task_id"],
                "purity_pct": cond["purity_pct"], "mixture_pct": cond["mixture_pct"],
                "region": name,
                "predicted_cn_bar": f"{pred:.4f}",
                "observed_cn_bar": f"{obs:.4f}",
                "abs_error": f"{obs - pred:.4f}",
                "pct_error": f"{100 * (obs - pred) / pred:.2f}",
            })

    if missing:
        print(f"\nWARNING: {len(missing)} condition(s) had no segmentation CSV:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ------------------------------------------------------------- summary
    print(f"\nwrote {args.out}  ({len(rows)} rows)")
    print("\nper-region error, over all conditions where the region applies:")
    print(f"{'region':<14} {'n':>4} {'mean %err':>10} {'min %err':>10} {'max %err':>10} {'max |abs|':>10}")
    by_region = defaultdict(list)
    for r in rows:
        by_region[r["region"]].append(r)
    for name in ["baseline", "A_clonal", "B_c2_private", "C_c3_private"]:
        rs = by_region.get(name)
        if not rs:
            continue
        pes = [float(r["pct_error"]) for r in rs]
        aes = [abs(float(r["abs_error"])) for r in rs]
        print(f"{name:<14} {len(rs):>4} {sum(pes)/len(pes):>10.2f} "
              f"{min(pes):>10.2f} {max(pes):>10.2f} {max(aes):>10.4f}")

    print("\ne2 vs e4 at matched conditions (e2 should sit at or below e4):")
    pair = defaultdict(dict)
    for r in rows:
        if r["experiment"] in ("e2", "e4") and r["region"] != "baseline":
            pair[(r["condition"], r["region"])][r["experiment"]] = float(r["observed_cn_bar"])
    both = [v for v in pair.values() if "e2" in v and "e4" in v]
    below = sum(1 for v in both if v["e2"] <= v["e4"])
    total = len(both)
    print(f"  e2 <= e4 in {below} of {total} matched (condition, region) comparisons")
    return 0


if __name__ == "__main__":
    sys.exit(main())
