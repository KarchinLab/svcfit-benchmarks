#!/usr/bin/env python3
"""Join SV -> clone assignment against the per-condition clone fractions.

Emits one row per (experiment, condition, SV) giving that SV's cellular
fraction: the fraction of cells in the simulated tumour that carry it.

On a hemizygous chrX a carrier cell has the SV on its only copy, so the
cellular fraction is exactly what SVCFit's corrected estimator should recover.
That is the whole point of the simulation.

  clonal      -> present in c2 and c3      -> fraction = purity
  c2_private  -> present in c2 only        -> fraction = c2 percentage
  c3_private  -> present in c3 only        -> fraction = c3 percentage

DELIBERATELY NOT EMITTED: copy number at the SV locus. cn_bar is estimated from
depth by 05_chrx_depth_segmentation_sim.R at scoring time; scoring is out of scope
here. Encoding an analytic copy number now would create a second, unvalidated
source of truth for the same quantity.

EXPERIMENT -> SV SET. e1 uses its own copy-neutral SV set. e2 and e4 share the
'cnv' set: the same SVs in the same reference coordinates, differing only in
whether the amplification was applied before or after them (01_chrx_hack.sh).
Cellular fraction does not depend on that ordering — clone membership does not
change — so both experiments take the same fractions.

COORDINATES are reference space for every set, verified by set-diffing the
assignment table against c2/c3/c22/c33.bed. e4's clone genomes were built from
aft_*.bed, whose coordinates are shifted by the upstream amplification, but the
assignment table records the reference-space position, which is what an aligned
BAM will report.
"""

import argparse
import os
import sys
from collections import defaultdict

# experiment -> which SV set in sv_subclone_assignment.tsv it was built from
EXPERIMENT_SV_SET = {"e1": "e1", "e2": "cnv", "e4": "cnv"}

# SV set -> (clonal bed, c2-carrying bed, c3-carrying bed), for the guard below
SET_BEDS = {
    "e1":  ("c1.bed",  "c2.bed",  "c3.bed"),
    "cnv": ("c11.bed", "c22.bed", "c33.bed"),
}

MEMBERSHIPS = ("clonal", "c2_private", "c3_private")


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def read_bed_coords(path):
    """(chrom, start, end) triples from a VISOR HACk bed."""
    coords = set()
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            coords.add((f[0], int(f[1]), int(f[2])))
    return coords


def fail(msg):
    print(f"  GUARD FAILED: {msg}", file=sys.stderr)
    return 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditions", required=True)
    ap.add_argument("--assignment", required=True)
    ap.add_argument("--sv-bed-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    conditions = read_tsv(args.conditions)
    svs = read_tsv(args.assignment)

    # ---------------------------------------------------------------- guards
    # The assignment table is the only record of which clone carries which SV.
    # If it disagreed with the beds that actually built the clone genomes, every
    # fraction below would be confidently wrong with nothing looking broken —
    # the same failure shape as a diploid normal.
    # So verify membership against the beds rather than trusting the table.
    print("=== guards ===")
    rc = 0

    by_set = defaultdict(list)
    for sv in svs:
        by_set[sv["set"]].append(sv)

    missing_sets = set(EXPERIMENT_SV_SET.values()) - set(by_set)
    if missing_sets:
        rc |= fail(f"assignment table lacks SV set(s): {sorted(missing_sets)}")

    for sv_set, (clonal_bed, c2_bed, c3_bed) in SET_BEDS.items():
        if sv_set not in by_set:
            continue
        rows = by_set[sv_set]

        bad = [r["membership"] for r in rows if r["membership"] not in MEMBERSHIPS]
        if bad:
            rc |= fail(f"set {sv_set}: unknown membership value(s) {sorted(set(bad))}")
            continue

        def coords(memberships):
            return {(r["CHROM"], int(r["start"]), int(r["end"]))
                    for r in rows if r["membership"] in memberships}

        checks = [
            (clonal_bed, coords({"clonal"}),                "clonal"),
            (c2_bed,     coords({"clonal", "c2_private"}),  "clonal+c2_private"),
            (c3_bed,     coords({"clonal", "c3_private"}),  "clonal+c3_private"),
        ]
        for bed_name, from_table, label in checks:
            bed_path = os.path.join(args.sv_bed_dir, bed_name)
            if not os.path.exists(bed_path):
                rc |= fail(f"missing {bed_path}")
                continue
            from_bed = read_bed_coords(bed_path)
            if from_table != from_bed:
                rc |= fail(
                    f"set {sv_set}: {label} ({len(from_table)} SVs) != {bed_name} "
                    f"({len(from_bed)} SVs); "
                    f"table-only={len(from_table - from_bed)} bed-only={len(from_bed - from_table)}")
            else:
                print(f"  {sv_set:>3} {label:<18} == {bed_name:<9} "
                      f"({len(from_bed)} SVs)  OK")

    # Clone fractions must partition the cell population exactly. decode_task()
    # uses integer arithmetic; a truncating division would silently lose cells.
    for c in conditions:
        total = int(c["c2_pct"]) + int(c["c3_pct"]) + int(c["normal_pct"])
        if total != 100:
            rc |= fail(f"{c['experiment']}/{c['condition']}: clone fractions "
                       f"sum to {total}%, not 100%")
        if int(c["c2_pct"]) + int(c["c3_pct"]) != int(c["purity_pct"]):
            rc |= fail(f"{c['experiment']}/{c['condition']}: c2+c3 != purity")
    print(f"  clone fractions sum to 100% in all {len(conditions)} conditions  OK")

    unknown_exp = {c["experiment"] for c in conditions} - set(EXPERIMENT_SV_SET)
    if unknown_exp:
        rc |= fail(f"condition grid has experiment(s) with no SV set mapping: "
                   f"{sorted(unknown_exp)}")

    if rc:
        print("\nRESULT: FAILED — ground truth not written", file=sys.stderr)
        return 1

    # ------------------------------------------------------------- stable ids
    # Ordinal within the SV set, by coordinate, so an SV keeps one id across
    # every condition and both experiments that share the set.
    sv_id = {}
    for sv_set, rows in by_set.items():
        ordered = sorted(rows, key=lambda r: (r["CHROM"], int(r["start"]), int(r["end"])))
        for i, r in enumerate(ordered, start=1):
            sv_id[(sv_set, r["CHROM"], r["start"], r["end"])] = f"{sv_set}_{i:03d}"

    # ------------------------------------------------------------------ join
    columns = ["experiment", "condition", "task_id", "purity_pct", "mixture_pct",
               "c2_pct", "c3_pct", "normal_pct",
               "sv_id", "CHROM", "start", "end", "sv_class", "sv_info",
               "membership", "carrier_clones",
               "cellular_fraction", "cellular_fraction_pct"]

    n = 0
    with open(args.out, "w") as out:
        out.write("\t".join(columns) + "\n")
        for cond in sorted(conditions, key=lambda c: int(c["task_id"])):
            sv_set = EXPERIMENT_SV_SET[cond["experiment"]]
            pct_for = {
                "clonal":     int(cond["purity_pct"]),
                "c2_private": int(cond["c2_pct"]),
                "c3_private": int(cond["c3_pct"]),
            }
            carriers_for = {"clonal": "c2,c3", "c2_private": "c2", "c3_private": "c3"}

            for sv in sorted(by_set[sv_set],
                             key=lambda r: (r["CHROM"], int(r["start"]), int(r["end"]))):
                pct = pct_for[sv["membership"]]
                out.write("\t".join([
                    cond["experiment"], cond["condition"], cond["task_id"],
                    cond["purity_pct"], cond["mixture_pct"],
                    cond["c2_pct"], cond["c3_pct"], cond["normal_pct"],
                    sv_id[(sv_set, sv["CHROM"], sv["start"], sv["end"])],
                    sv["CHROM"], sv["start"], sv["end"],
                    sv["class"], sv["info"],
                    sv["membership"], carriers_for[sv["membership"]],
                    f"{pct / 100:.4f}", str(pct),
                ]) + "\n")
                n += 1

    print(f"\nwrote {args.out}")
    print(f"  {n} rows = {len(conditions)} conditions x SVs per experiment")
    print("RESULT: ground truth written and verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
