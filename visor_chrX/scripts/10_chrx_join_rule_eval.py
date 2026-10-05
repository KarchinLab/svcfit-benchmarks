#!/usr/bin/env python3
"""Score the three segment-to-SV join rules against ANALYTIC truth.

The three rules, all implemented in 08_chrx_segment_sv_join.py:

    unweighted   mean of the two breakpoint segments' cn_bar   (the recorded ruling)
    weighted     length-weighted mean over the SV span         (--weighted)
    majority     single segment of maximum overlap, wholesale  (--majority)

The third is the AUTOSOMAL rule. SVCFit joins FACETS segments to SVs with
`slice_max(width, with_ties = FALSE)` in assign_cnv.R and again in
assign_bkg_cnv.R.

CAVEAT ON THE TRUTH USED HERE, stated because 04_chrx_ground_truth.py:15
deliberately refused to emit it: "DELIBERATELY NOT EMITTED: copy number at the
SV locus ... Encoding an analytic copy number now would create a second,
unvalidated truth." That refusal was right for the ground-truth table, which
feeds scoring. This script needs a per-SV copy number to score the JOIN, so it
builds one -- and therefore validates it against the two independently verified
anchors in the docs before reporting anything (--validate-only shows them).

THE MODEL. cn_bar is the mean copies of the locus per cell, over the SV's own
span, which is the definition the depth-based cn_bar estimator uses through its recut
worked example. Per clone k with amplification factor a(p) from dup_c{22,33}.bed
(reference coordinates), for a position p inside an SV the clone carries:

  e1  no amplification, a = 1
      deletion  -> 0        inversion -> 1        tandem dup (r=2) -> 2

  e4  AMPLIFICATION FIRST, then SV  (01_chrx_hack.sh:131-137)
      the SV acts on one copy of an already-amplified region
      deletion  -> a - 1    inversion -> a

  e2  SV FIRST, then amplification  (01_chrx_hack.sh:118-128)
      the SV acts on a single-copy genome; the amp then multiplies what remains
      deletion  -> 0        inversion -> a

Then cn_bar(p) = f_c2*cn_c2(p) + f_c3*cn_c3(p) + f_n*1, and the SV's true
cn_bar is the mean of that over its span (the amp edge can fall inside a span,
so this integrates piecewise rather than sampling the midpoint).

The e2/e4 asymmetry is why the segmentation eval checks e2 <= e4: at an SV locus
inside an amplification, SV-first drives copies to 0 while amp-first leaves a-1.
No SV overlaps another in any condition (verified), so each span is affected
only by its own SV and the amplifications.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

CHRX_LEN = 156_040_895
CLONES = ("c2", "c3")


def read_dup_bed(path):
    """{(start, end): copies} from a VISOR tandem-duplication bed."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 5:
                continue
            out[(int(f[1]), int(f[2]))] = int(f[4])
    return out


def amp_breaks(amps):
    b = set()
    for (s, e) in amps:
        b.add(s)
        b.add(e)
    return b


def amp_at(amps, pos):
    for (s, e), cn in amps.items():
        if s <= pos < e:
            return cn
    return 1


def clone_cn(exp, a, carries, sv_class, r):
    """Copies of a position in one clone, given its amplification factor."""
    if not carries:
        return a
    if sv_class == "inversion":
        return a
    if sv_class == "deletion":
        return 0 if exp == "e2" else a - 1
    if sv_class == "tandem duplication":
        return a * r if exp == "e2" else a + (r - 1)
    raise ValueError(f"unhandled sv_class {sv_class!r}")


def true_cn_bar(exp, start, end, carriers, sv_class, r, f2, f3, fn, amp):
    """Mean copies per cell over [start, end), integrated piecewise."""
    cuts = sorted({start, end} | {b for b in amp["breaks"] if start < b < end})
    num = 0.0
    for lo, hi in zip(cuts, cuts[1:]):
        mid = (lo + hi) // 2
        cn = fn * 1.0
        for clone, f in (("c2", f2), ("c3", f3)):
            a = 1 if exp == "e1" else amp[clone].get("at", lambda p: 1)(mid)
            cn += f * clone_cn(exp, a, clone in carriers, sv_class, r)
        num += cn * (hi - lo)
    return num / (end - start)


def load_joined(path):
    """{(start): cn_bar} from a Stage 4 <sample>_cn_bar.csv."""
    out = {}
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[int(r["POS"])] = None if r["cn_bar"] == "NA" else float(r["cn_bar"])
    return out


def summarize(name, diffs):
    if not diffs:
        return f"{name:<12} {'no rows':>10}"
    n = len(diffs)
    absd = sorted(abs(d) for d in diffs)
    mean = sum(absd) / n
    med = absd[n // 2]
    p95 = absd[int(0.95 * (n - 1))]
    within = sum(1 for d in absd if d <= 0.05)
    return (f"{name:<12} {n:>6} {mean:>10.4f} {med:>10.4f} {p95:>10.4f} "
            f"{absd[-1]:>10.4f} {within:>7} ({100*within/n:>5.1f}%)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--truth", default=None)
    ap.add_argument("--sv-bed-dir", default=None)
    ap.add_argument("--rules", default="unweighted=cn_bar10k,weighted=cn_bar10k_w,majority=cn_bar10k_m",
                    help="comma-separated name=dir")
    ap.add_argument("--out", default=None, help="per-SV TSV of truth and all rules")
    ap.add_argument("--validate-only", action="store_true")
    args = ap.parse_args()

    base = args.base
    # Coverage-named; see the note in chrx_common.sh. REQUIRED, no default.
    cov = os.environ.get("COV", "")
    if not cov:
        raise SystemExit("COV is not set, and there is no default. Run e.g. COV=50 ...")
    ctag = "c" + cov
    truth_path = args.truth or os.path.join(
        base, "truth", f"chrx_sv_ground_truth_{ctag}.tsv")
    bed_dir = args.sv_bed_dir or os.path.join(base, "truth", "sv_beds")

    amp = {}
    for clone, bed in (("c2", "dup_c22.bed"), ("c3", "dup_c33.bed")):
        d = read_dup_bed(os.path.join(bed_dir, bed))
        amp[clone] = {"d": d, "at": (lambda dd: (lambda p: amp_at(dd, p)))(d)}
    amp["breaks"] = amp_breaks(amp["c2"]["d"]) | amp_breaks(amp["c3"]["d"])

    print("amplification regions (reference coordinates):")
    for clone in CLONES:
        for (s, e), cn in sorted(amp[clone]["d"].items()):
            print(f"  {clone}  {s:>11,}-{e:>11,}  copies={cn}")

    rows = list(csv.DictReader(open(truth_path), delimiter="\t"))
    by_cond = defaultdict(list)
    for r in rows:
        by_cond[(r["experiment"], r["condition"])].append(r)

    # ---------------- validate the model against the verified anchors --------
    print("\nvalidation against anchors verified elsewhere in the docs:")
    # The condition carries the coverage tag, but the quantity being checked does
    # not: cn_bar here is a copy number built from purity and mixture, and those
    # are coverage-independent by construction, so the anchors do not name a coverage.
    anchors = [
        # (exp, cond, sv start, expected, source)
        ("e4", f"{ctag}p80m10", 7534655, 1.80,
         "segmenter read cn_bar 1.8103"),
    ]
    ok = True
    for exp, cond, start, expect, src in anchors:
        hit = [r for r in by_cond[(exp, cond)] if int(r["start"]) == start]
        if not hit:
            print(f"  MISS  {exp}/{cond} start={start:,} not in truth table")
            ok = False
            continue
        r = hit[0]
        got = true_cn_bar(exp, int(r["start"]), int(r["end"]),
                          set(r["carrier_clones"].split(",")), r["sv_class"],
                          int(r["sv_info"]) if r["sv_info"] not in ("None", "") else 1,
                          int(r["c2_pct"]) / 100, int(r["c3_pct"]) / 100,
                          int(r["normal_pct"]) / 100, amp)
        flag = "OK  " if abs(got - expect) < 1e-6 else "FAIL"
        ok &= abs(got - expect) < 1e-6
        print(f"  {flag}  {exp}/{cond} {r['sv_class']:<18} model={got:.4f} "
              f"expected={expect:.4f}   [{src}]")

    # e2 <= e4 at SV loci, the asymmetry the ordering implies
    pairs = viol = 0
    for (exp, cond), rs in by_cond.items():
        if exp != "e2":
            continue
        e4 = {int(x["start"]): x for x in by_cond[("e4", cond)]}
        for r in rs:
            s = int(r["start"])
            if s not in e4:
                continue
            kw = dict(f2=int(r["c2_pct"]) / 100, f3=int(r["c3_pct"]) / 100,
                      fn=int(r["normal_pct"]) / 100, amp=amp)
            a = true_cn_bar("e2", s, int(r["end"]), set(r["carrier_clones"].split(",")),
                            r["sv_class"], 1, **kw)
            b = true_cn_bar("e4", s, int(e4[s]["end"]), set(e4[s]["carrier_clones"].split(",")),
                            e4[s]["sv_class"], 1, **kw)
            pairs += 1
            viol += (a > b + 1e-9)
    print(f"  {'OK  ' if not viol else 'FAIL'}  e2 <= e4 at matched SV loci: "
          f"{pairs - viol}/{pairs}")
    ok &= not viol

    if not ok:
        print("\nMODEL FAILED VALIDATION -- not scoring.", file=sys.stderr)
        return 1
    if args.validate_only:
        return 0

    # ---------------- score each rule ---------------------------------------
    rules = [kv.split("=", 1) for kv in args.rules.split(",")]
    joined = {}
    for name, d in rules:
        joined[name] = {k: load_joined(os.path.join(base, d, f"{k[0]}_{k[1]}_cn_bar.csv"))
                        for k in by_cond}

    diffs = {name: [] for name, _ in rules}
    diffs_by = {name: defaultdict(list) for name, _ in rules}
    out_rows = []
    missing = defaultdict(int)

    for (exp, cond), rs in sorted(by_cond.items()):
        for r in rs:
            s, e = int(r["start"]), int(r["end"])
            t = true_cn_bar(exp, s, e, set(r["carrier_clones"].split(",")), r["sv_class"],
                            int(r["sv_info"]) if r["sv_info"] not in ("None", "") else 1,
                            int(r["c2_pct"]) / 100, int(r["c3_pct"]) / 100,
                            int(r["normal_pct"]) / 100, amp)
            row = {"experiment": exp, "condition": cond, "sv_id": r["sv_id"],
                   "CHROM": r["CHROM"], "start": s, "end": e,
                   "sv_class": r["sv_class"], "membership": r["membership"],
                   "true_cn_bar": f"{t:.4f}"}
            for name, _ in rules:
                j = joined[name][(exp, cond)]
                v = j.get(s) if j else None
                row[name] = "NA" if v is None else f"{v:.4f}"
                if v is None:
                    missing[name] += 1
                    continue
                d = v - t
                diffs[name].append(d)
                diffs_by[name][exp].append(d)
                diffs_by[name][r["sv_class"]].append(d)
            out_rows.append(row)

    hdr = (f"\n{'rule':<12} {'n':>6} {'mean|err|':>10} {'median':>10} {'p95':>10} "
           f"{'max':>10} {'within 0.05':>15}")
    print("\n=== all SVs, all 45 conditions ===")
    print(hdr)
    for name, _ in rules:
        print(summarize(name, diffs[name]))
        if missing[name]:
            print(f"{'':<12} ({missing[name]} SVs absent from {name} output)")

    for split in ("e1", "e2", "e4", "deletion", "inversion", "tandem duplication"):
        if not any(diffs_by[n][split] for n, _ in rules):
            continue
        print(f"\n=== {split} ===")
        print(hdr)
        for name, _ in rules:
            print(summarize(name, diffs_by[name][split]))

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out_rows[0].keys()))
            w.writeheader()
            w.writerows(out_rows)
        print(f"\nwrote {args.out} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
