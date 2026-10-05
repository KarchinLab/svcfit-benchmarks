#!/usr/bin/env python3
"""Segment-to-SV join: per-segment cn_bar -> per-SV cn_bar.

    RULE: cn_bar for an SV is the cn_bar of the single segment of
    MAXIMUM OVERLAP with the SV span, taken wholesale.

This is the rule SVCFit uses on the autosomes: assign_cnv.R does
slice_max(width, with_ties = FALSE) over the SV span, and assign_bkg_cnv.R does
the same per flank.

Scoring all three rules against analytic truth (10_chrx_join_rule_eval.py) puts this one
first on every stratum, so consistency and accuracy agree here:

    rule         mean |err|   within 0.05
    unweighted      0.1535        52.3%     <- the superseded ruling
    weighted        0.0391        74.7%
    majority        0.0331        80.7%     <- default

Nothing lands in hemizygous_cnv_unresolved on account of a boundary: every
spanning SV is still estimated, as the superseded ruling also intended.
--weighted and --unweighted keep the alternatives runnable for comparison.

OUTPUT is the Stage 4 contract exactly: <sample>_cn_bar.csv with columns
CHROM, POS, cn_bar, one row per hemizygous SV. Diagnostics go to a SEPARATE
audit file so the primary output stays strictly to spec for whatever reads it.

SEGMENT COVERAGE. DNAcopy reports loc.start and loc.end as BIN START positions,
so a segment's loc.end is the start of its last bin, not the end of it. Taken
literally the segments leave 100 kb holes between them -- and real breakpoints
land in those holes: chrX:7,533,655 sits between a segment ending 7,500,001 and
the next starting 7,600,001. Each segment is therefore extended to the next
segment's start (and the last to the chromosome end), which is identical to
adding one bin width to loc.end. Without this, ~every amplification-boundary SV
would report no_segment.
"""

import argparse
import csv
import os
import sys
from collections import Counter

# hg38 chrX. ASSEMBLY-DEPENDENT: the prostate mixture benchmark is GRCh37/hs37d5, whose X is
# 155,270,560 -- 770,335 bp shorter. Overridable with --chrom-len; the default keeps the COMBAT
# cohort and the chrX simulation unchanged.
#
# It is used in one place, to extend the LAST segment to the end of the chromosome so a breakpoint
# past the final segment boundary still lands somewhere. Getting it too LONG, which is what the
# hg38 default does on GRCh37, is harmless -- the extra span holds no SVs. Getting it too SHORT
# would silently drop breakpoints beyond it, since find_segment() returns None and the SV is
# recorded as having no depth rather than as an error. Parameterised so neither has to be reasoned
# about again.
CHRX_LEN_HG38 = 156_040_895


def read_segments(path, chrom_len):
    """Sorted, gap-filled segments: [(start, end_exclusive, cn_bar), ...]."""
    segs = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            segs.append([int(r["start"]), int(r["end"]), float(r["cn_bar"])])
    if not segs:
        return []
    segs.sort(key=lambda s: s[0])
    out = []
    for i, (s, e, cn) in enumerate(segs):
        end = segs[i + 1][0] if i + 1 < len(segs) else max(e, chrom_len)
        out.append((s, end, cn))
    return out


def find_segment(segs, pos):
    for s, e, cn in segs:
        if s <= pos < e:
            return (s, e, cn)
    return None


def read_svs(path, chrom_col, start_col, end_col, chrom_filter):
    with open(path) as fh:
        sample = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in sample else ","
        rows = []
        for r in csv.DictReader(fh, delimiter=delim):
            if chrom_filter and r[chrom_col] != chrom_filter:
                continue
            rows.append((r[chrom_col], int(r[start_col]), int(r[end_col]), r))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--segments", required=True, help="chrx_seg_*.csv from script 05")
    ap.add_argument("--svs", required=True, help="TSV/CSV with chrom/start/end columns")
    ap.add_argument("--out", required=True, help="<sample>_cn_bar.csv (CHROM,POS,cn_bar)")
    ap.add_argument("--audit", help="optional per-SV diagnostics TSV")
    ap.add_argument("--bg-out", dest="bg_out",
                    help="optional <sample>_bg_cn.csv (CHROM,POS,bg_cn): the flanking copy number, "
                         "which hemizygous_del_svcf() needs as kappa for the CNV-first form")
    ap.add_argument("--chrom-col", default="CHROM")
    ap.add_argument("--start-col", default="start")
    ap.add_argument("--end-col", default="end")
    ap.add_argument("--chrom", default="chrX", help="restrict to this contig")
    ap.add_argument("--out-chrom", dest="out_chrom", default=None,
                    help="write this contig name in the output instead of the input's. Needed when "
                         "the consumer normalises names: SVCFit's load_data() rewrites every CHROM "
                         "to paste0('chr', sub('^chr','',CHROM)), so a GRCh37 run whose VCFs say 'X' "
                         "is matched internally as 'chrX'. calc_svcf looks up hemi_cn_bar on "
                         "(CHROM, POS), so a table saying 'X' would miss every row and every SV "
                         "would silently fall through as having no depth. Default: unchanged.")
    ap.add_argument("--chrom-len", dest="chrom_len", type=int, default=CHRX_LEN_HG38,
                    help="length of that contig in bp, used only to extend the last segment to "
                         "the chromosome end (default %(default)d, hg38 chrX; GRCh37 X is 155270560)")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--majority", action="store_true",
                      help="DEFAULT. Single segment of maximum overlap with the "
                           "SV span, taken wholesale -- the autosomal rule "
                           "(assign_cnv.R:41-44)")
    mode.add_argument("--weighted", action="store_true",
                      help="length-weighted mean of segment cn_bar over the SV span")
    mode.add_argument("--unweighted", action="store_true",
                      help="unweighted mean of the two breakpoint segments' cn_bar "
                           "-- kept for comparison")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    segs = read_segments(args.segments, args.chrom_len)
    if not segs:
        print(f"ERROR: no segments in {args.segments}", file=sys.stderr)
        return 1
    svs = read_svs(args.svs, args.chrom_col, args.start_col, args.end_col, args.chrom)
    # An empty SV set after filtering is almost always a contig-name mismatch -- SVs on "X" against
    # --chrom chrX, or the reverse -- and it would otherwise produce a valid, empty output file.
    if not svs:
        print(f"ERROR: no SVs on contig {args.chrom!r} in {args.svs}. Check --chrom against the "
              f"{args.chrom_col} column; GRCh37 uses 'X' where hg38 uses 'chrX'.", file=sys.stderr)
        return 1

    results, status = [], Counter()
    for chrom, start, end, _raw in svs:
        s1 = find_segment(segs, start)
        s2 = find_segment(segs, end)

        if s1 is None or s2 is None:
            # Never silently estimated. A breakpoint outside every segment means
            # the segmentation had no informative bin there.
            cn, st, detail = None, "no_segment", ""
        elif s1 == s2:
            cn, st = s1[2], "single_segment"
            detail = f"{s1[2]:.4f}"
        elif args.weighted:
            num = den = 0.0
            for a, b, c in segs:
                ov = min(b, end) - max(a, start)
                if ov > 0:
                    num += c * ov
                    den += ov
            cn, st = num / den, "spanning_weighted"
            detail = f"{s1[2]:.4f}|{s2[2]:.4f}"
        elif args.unweighted:
            cn, st = (s1[2] + s2[2]) / 2, "spanning_mean"
            detail = f"{s1[2]:.4f}|{s2[2]:.4f}"
        else:
            # assign_cnv.R:41-44 -- slice_max(width, with_ties = FALSE) over the
            # SV span. Ties go to the leftmost segment, which is what with_ties
            # = FALSE does given findOverlaps' ascending subjectHits order.
            best_ov, cn = 0, None
            for a, b, c in segs:
                ov = min(b, end) - max(a, start)
                if ov > best_ov:
                    best_ov, cn = ov, c
            st = "spanning_majority"
            detail = f"{s1[2]:.4f}|{s2[2]:.4f}"

        status[st] += 1
        results.append({
            # Both outputs and the audit take this, so --out-chrom cannot rename one and not another.
            "CHROM": (args.out_chrom or chrom), "start": start, "end": end,
            "cn_bar": cn, "status": st, "breakpoint_cn_bars": detail,
            "n_segments_spanned": sum(1 for a, b, _ in segs if min(b, end) > max(a, start)),
        })

    # ---- background copy number, for the CNV-first deletion form ----
    #
    # kappa: the copies the locus would have WITHOUT the deletion. hemizygous_del_svcf() needs it,
    # because 1 - cn_bar assumes kappa = 1 and is biased by -(kappa - 1) when it is not
    # (Supplementary Note S5).
    #
    # DEFINITION IS 13_chrx_score_svcfit.R's attach_background_cn(), transcribed rather than
    # reinvented: the length-weighted mean cn_bar over the SV span extended by FLANK_BP on each
    # side, MINUS the SV's own span. Deriving it a second way here would mean the cohort and the
    # simulation that validated the form were measuring different quantities.
    FLANK_BP = 500_000

    def background_cn(start, end):
        lo, hi = start - FLANK_BP, end + FLANK_BP
        num = den = 0.0
        for s, e, cn in segs:
            flank_w = max(0, min(e, hi) - max(s, lo)) - max(0, min(e, end) - max(s, start))
            if flank_w > 0:
                num += cn * flank_w
                den += flank_w
        return num / den if den > 0 else None

    for r in results:
        r["bg_cn"] = background_cn(int(r["start"]), int(r["end"]))

    if args.bg_out:
        with open(args.bg_out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["CHROM", "POS", "bg_cn"])
            for r in results:
                w.writerow([r["CHROM"], r["start"],
                            "NA" if r["bg_cn"] is None else f"{r['bg_cn']:.4f}"])

    # ---- primary output: the Stage 4 contract, exactly three columns ----
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["CHROM", "POS", "cn_bar"])
        for r in results:
            w.writerow([r["CHROM"], r["start"],
                        "NA" if r["cn_bar"] is None else f"{r['cn_bar']:.4f}"])

    if args.audit:
        with open(args.audit, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(results[0].keys()))
            w.writeheader()
            for r in results:
                r = dict(r)
                r["cn_bar"] = "NA" if r["cn_bar"] is None else f"{r['cn_bar']:.4f}"
                r["bg_cn"]  = "NA" if r.get("bg_cn") is None else f"{r['bg_cn']:.4f}"
                w.writerow(r)

    if not args.quiet:
        print(f"{os.path.basename(args.out)}: {len(results)} SVs, "
              f"{len(segs)} segments — " +
              ", ".join(f"{k}={v}" for k, v in sorted(status.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
