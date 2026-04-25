#!/usr/bin/env python3
"""Compare sorted BED3 peak/region sets (MACS3 bdgpeakcall vs C++ diagnostic)."""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass(frozen=True)
class Iv:
    chrom: str
    start: int
    end: int


def read_bed3(path: str) -> List[Iv]:
    rows: List[Iv] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            c, s, e = parts[0], int(parts[1]), int(parts[2])
            if e > s:
                rows.append(Iv(c, s, e))
    rows.sort(key=lambda x: (x.chrom, x.start, x.end))
    return rows


def stats(ivs: List[Iv]) -> Tuple[int, int, float]:
    if not ivs:
        return 0, 0, 0.0
    widths = [x.end - x.start for x in ivs]
    tot = sum(widths)
    ws = sorted(widths)
    m = len(ws) // 2
    if len(ws) % 2:
        med = float(ws[m])
    else:
        med = (ws[m - 1] + ws[m]) / 2.0
    return len(ivs), tot, med


def jaccard(a: List[Iv], b: List[Iv]) -> float:
    """Unweighted interval Jaccard on union of chromosomes (merge overlaps within each set)."""
    def merged_len(ivs: List[Iv]) -> int:
        if not ivs:
            return 0
        tot = 0
        cur_s, cur_e = ivs[0].start, ivs[0].end
        cur_c = ivs[0].chrom
        for x in ivs[1:]:
            if x.chrom == cur_c and x.start < cur_e:
                cur_e = max(cur_e, x.end)
            else:
                tot += cur_e - cur_s
                cur_c, cur_s, cur_e = x.chrom, x.start, x.end
        tot += cur_e - cur_s
        return tot

    def inter_len(aa: List[Iv], bb: List[Iv]) -> int:
        i = j = 0
        tot = 0
        while i < len(aa) and j < len(bb):
            x, y = aa[i], bb[j]
            if x.chrom < y.chrom:
                i += 1
            elif y.chrom < x.chrom:
                j += 1
            else:
                lo = max(x.start, y.start)
                hi = min(x.end, y.end)
                if hi > lo:
                    tot += hi - lo
                if x.end <= y.end:
                    i += 1
                else:
                    j += 1
        return tot

    la = merged_len(a)
    lb = merged_len(b)
    inter = inter_len(a, b)
    union = la + lb - inter
    if union <= 0:
        return 1.0 if not a and not b else 0.0
    return inter / union


def reciprocal_any_overlap(a: List[Iv], b: List[Iv]) -> Tuple[float, float]:
    def count_any(xlist: List[Iv], ylist: List[Iv]) -> int:
        if not xlist:
            return 0
        j = 0
        n = 0
        for x in xlist:
            while j < len(ylist) and (ylist[j].chrom < x.chrom or (
                    ylist[j].chrom == x.chrom and ylist[j].end <= x.start)):
                j += 1
            k = j
            hit = False
            while k < len(ylist) and ylist[k].chrom == x.chrom and ylist[k].start < x.end:
                if min(x.end, ylist[k].end) > max(x.start, ylist[k].start):
                    hit = True
                    break
                k += 1
            if hit:
                n += 1
        return n

    ca = count_any(a, b)
    cb = count_any(b, a)
    fa = ca / len(a) if a else 1.0
    fb = cb / len(b) if b else 1.0
    return fa, fb


def greedy_match_dists(a: List[Iv], b: List[Iv]) -> Tuple[List[int], List[int], int, int]:
    """Greedy: for each interval in A, pick first B on same chrom with any overlap; record |dstart|, |dend|."""
    ds: List[int] = []
    de: List[int] = []
    used_b = set()
    j0 = 0
    for i, x in enumerate(a):
        j = j0
        best: Optional[Tuple[int, int, int]] = None  # idx, overlap, dist penalty
        while j < len(b):
            y = b[j]
            if y.chrom < x.chrom:
                j += 1
                j0 = j
                continue
            if y.chrom > x.chrom:
                break
            if y.end <= x.start:
                j += 1
                continue
            if y.start >= x.end:
                break
            ov = min(x.end, y.end) - max(x.start, y.start)
            if j not in used_b:
                pen = abs(x.start - y.start) + abs(x.end - y.end)
                if best is None or ov > best[1] or (ov == best[1] and pen < best[2]):
                    best = (j, ov, pen)
            j += 1
        if best is not None:
            idx = best[0]
            used_b.add(idx)
            y = b[idx]
            ds.append(abs(x.start - y.start))
            de.append(abs(x.end - y.end))
    matched = len(ds)
    only_a = len(a) - matched
    only_b = len(b) - len(used_b)
    return ds, de, only_a, only_b


def summarize_ints(xs: List[int]) -> str:
    if not xs:
        return "na"
    xs = sorted(xs)
    m = len(xs) // 2
    med = xs[m] if len(xs) % 2 else (xs[m - 1] + xs[m]) / 2.0
    return f"min={xs[0]} median={med} max={xs[-1]} mean={sum(xs)/len(xs):.4g}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("bed_a", help="BED3 or narrowPeak (first 3 cols used if more)")
    ap.add_argument("bed_b")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    args = ap.parse_args()
    a = read_bed3(args.bed_a)
    b = read_bed3(args.bed_b)
    lines = []
    ident = a == b
    lines.append(f"bed3_sorted_identical\t{ident}")
    na, ta, ma = stats(a)
    nb, tb, mb = stats(b)
    lines.append(f"n_regions_{args.label_a}\t{na}")
    lines.append(f"n_regions_{args.label_b}\t{nb}")
    lines.append(f"total_bp_{args.label_a}\t{ta}")
    lines.append(f"total_bp_{args.label_b}\t{tb}")
    lines.append(f"median_width_{args.label_a}\t{ma}")
    lines.append(f"median_width_{args.label_b}\t{mb}")
    jac = jaccard(a, b)
    lines.append(f"jaccard_interval_bp\t{jac:.12g}")
    ra, rb = reciprocal_any_overlap(a, b)
    lines.append(f"reciprocal_any_overlap_fraction_{args.label_a}\t{ra:.12g}")
    lines.append(f"reciprocal_any_overlap_fraction_{args.label_b}\t{rb:.12g}")
    ds, de, ua, ub = greedy_match_dists(a, b)
    lines.append(f"n_unique_to_{args.label_a}\t{ua}")
    lines.append(f"n_unique_to_{args.label_b}\t{ub}")
    lines.append(f"boundary_start_abs_diff_{args.label_a}_vs_{args.label_b}\t{summarize_ints(ds)}")
    lines.append(f"boundary_end_abs_diff_{args.label_a}_vs_{args.label_b}\t{summarize_ints(de)}")
    out = "\n".join(lines) + "\n"
    sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
