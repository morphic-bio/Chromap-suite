#!/usr/bin/env python3
"""Compare MACS3 callpeak narrowPeak/summits vs C++ diagnostic (Phase 6 metrics)."""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass
class Np:
    chrom: str
    start: int
    end: int
    name: str
    score: int
    strand: str
    signal: float
    pval: float
    qval: float
    peak: int
    line_no: int


@dataclass
class Smt:
    chrom: str
    start: int
    end: int
    name: str
    sc: float
    line_no: int


def read_narrow(path: str) -> List[Np]:
    r: List[Np] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for k, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browse"):
                continue
            p = line.split("\t")
            if len(p) < 10:
                continue
            r.append(
                Np(
                    p[0],
                    int(p[1]),
                    int(p[2]),
                    p[3],
                    int(float(p[4])),
                    p[5],
                    float(p[6]),
                    float(p[7]),
                    float(p[8]),
                    int(p[9]),
                    k,
                )
            )
    r.sort(key=lambda x: (x.chrom, x.start, x.end, x.peak))
    return r


def read_smt(path: str) -> List[Smt]:
    r: List[Smt] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for k, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            p = line.split("\t")
            if len(p) < 3:
                continue
            r.append(
                Smt(
                    p[0],
                    int(p[1]),
                    int(p[2]),
                    p[3] if len(p) > 3 else ".",
                    float(p[4]) if len(p) > 4 else 0.0,
                    k,
                )
            )
    r.sort(key=lambda x: (x.chrom, x.start, x.end))
    return r


def bed3_key(x: Np) -> Tuple[str, int, int]:
    return (x.chrom, x.start, x.end)


def stats_iv(xs: List[Tuple[str, int, int]]) -> Tuple[int, int, float]:
    if not xs:
        return 0, 0, 0.0
    w = [a[2] - a[1] for a in xs]
    tot = sum(w)
    ws = sorted(w)
    m = len(ws) // 2
    if len(ws) % 2:
        med = float(ws[m])
    else:
        med = (ws[m - 1] + ws[m]) / 2.0
    return len(xs), tot, med


def jaccard(a: List[Tuple[str, int, int]], b: List[Tuple[str, int, int]]) -> float:
    def inter_len(aa, bb) -> int:
        i = j = 0
        tot = 0
        while i < len(aa) and j < len(bb):
            x, y = aa[i], bb[j]
            if x[0] < y[0]:
                i += 1
            elif y[0] < x[0]:
                j += 1
            else:
                lo = max(x[1], y[1])
                hi = min(x[2], y[2])
                if hi > lo:
                    tot += hi - lo
                if x[2] <= y[2]:
                    i += 1
                else:
                    j += 1
        return tot

    def merge_len(ss: List[Tuple[str, int, int]]) -> int:
        if not ss:
            return 0
        tot = 0
        c, s, e = ss[0][0], ss[0][1], ss[0][2]
        for t in ss[1:]:
            if t[0] == c and t[1] < e:
                e = max(e, t[2])
            else:
                tot += e - s
                c, s, e = t[0], t[1], t[2]
        tot += e - s
        return tot

    a = sorted(a)
    b = sorted(b)
    la, lb = merge_len(a), merge_len(b)
    it = inter_len(a, b)
    u = la + lb - it
    if u <= 0:
        return 1.0 if not a and not b else 0.0
    return it / u


def pearson_r(xs: List[float], ys: List[float]) -> float:
    n = min(len(xs), len(ys))
    if n < 2:
        return float("nan")
    xs = xs[:n]
    ys = ys[:n]
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    d1 = sum((x - mx) ** 2 for x in xs) ** 0.5
    d2 = sum((y - my) ** 2 for y in ys) ** 0.5
    if d1 < 1e-30 or d2 < 1e-30:
        return float("nan")
    return num / (d1 * d2)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--macs3-np", required=True, help="MACS3 *_peaks.narrowPeak")
    ap.add_argument("--cpp-np", required=True, help="C++ diagnostic narrowPeak")
    ap.add_argument("--macs3-summits", default="", help="Optional MACS3 *_summits.bed")
    ap.add_argument("--cpp-summits", default="", help="Optional C++ summits BED5")
    ap.add_argument("--out-tsv", required=True, help="Summary TSV")
    args = ap.parse_args()

    m = read_narrow(args.macs3_np)
    c = read_narrow(args.cpp_np)
    a_iv = [bed3_key(x) for x in m]
    b_iv = [bed3_key(x) for x in c]
    ident = a_iv == b_iv
    n_a, t_a, w_a = stats_iv(a_iv)
    n_b, t_b, w_b = stats_iv(b_iv)
    jac = jaccard(a_iv, b_iv)

    lines: List[str] = []
    lines.append("metric\tvalue")
    lines.append(f"bed3_sorted_identical\t{1 if ident else 0}")
    lines.append(f"n_intervals_macs3\t{n_a}")
    lines.append(f"n_intervals_cpp\t{n_b}")
    lines.append(f"total_bp_macs3\t{t_a}")
    lines.append(f"total_bp_cpp\t{t_b}")
    lines.append(f"median_width_macs3\t{w_a}")
    lines.append(f"median_width_cpp\t{w_b}")
    lines.append(f"jaccard_bed3_intervals\t{jac:.12g}")

    # match intervals (greedy: exact key first + overlap for rest)
    used_c = set()
    pairs: List[Tuple[Np, Np]] = []
    c_key_to_i = {bed3_key(x): i for i, x in enumerate(c)}
    matched_mi = set()
    for mi, x in enumerate(m):
        bk = bed3_key(x)
        if bk in c_key_to_i and c_key_to_i[bk] not in used_c:
            j = c_key_to_i[bk]
            used_c.add(j)
            matched_mi.add(mi)
            pairs.append((x, c[j]))
    for mi, x in enumerate(m):
        if mi in matched_mi:
            continue
        j = 0
        best: Optional[Tuple[int, int, int]] = None
        while j < len(c):
            y = c[j]
            if y.chrom < x.chrom:
                j += 1
            elif y.chrom > x.chrom:
                break
            if y.end <= x.start:
                j += 1
            elif y.start >= x.end:
                break
            else:
                if j not in used_c:
                    ov = min(x.end, y.end) - max(x.start, y.start)
                    pen = abs(x.start - y.start) + abs(x.end - y.end)
                    if best is None or ov > best[0] or (ov == best[0] and pen < best[1]):
                        best = (ov, pen, j)
            j += 1
        if best is not None:
            j = best[2]
            used_c.add(j)
            matched_mi.add(mi)
            pairs.append((x, c[j]))
    n_match = len(pairs)
    only_m = [m[i] for i in range(len(m)) if i not in matched_mi]
    only_c = [c[i] for i in range(len(c)) if i not in used_c]
    lines.append(f"n_intervals_matched_greedy\t{n_match}")
    lines.append(f"n_only_macs3\t{len(only_m)}")
    lines.append(f"n_only_cpp\t{len(only_c)}")

    ra = n_match / len(m) if m else 0.0
    rb = n_match / len(c) if c else 0.0
    lines.append(f"reciprocal_any_overlap_fraction_macs3\t{ra:.12g}")
    lines.append(f"reciprocal_any_overlap_fraction_cpp\t{rb:.12g}")

    off_eq = 0
    sdist: List[int] = []
    s_sig: List[Tuple[float, float, float, float, float, float]] = []
    for x, y in pairs:
        if x.peak == y.peak:
            off_eq += 1
        absm = (x.start + x.peak) if x.start + x.peak >= 0 else 0
        absc = (y.start + y.peak) if y.start + y.peak >= 0 else 0
        sdist.append(abs(absm - absc))
        s_sig.append((x.signal, y.signal, x.pval, y.pval, x.qval, y.qval))
    lines.append(
        f"peak_offset_equal_count\t{off_eq} / {n_match}" if n_match else "peak_offset_equal_count\t0 / 0"
    )
    if sdist:
        sdist.sort()
        k = len(sdist) // 2
        med = sdist[k] if len(sdist) % 2 else (sdist[k - 1] + sdist[k]) / 2.0
        lines.append(f"summit_1bp_distance_abs_median_bp\t{med}")
        lines.append(
            f"summit_1bp_distance_abs_min_median_max\t{min(sdist)}\t{med}\t{max(sdist)}"
        )
    else:
        lines.append("summit_1bp_distance_abs_median_bp\tna")

    if s_sig and n_match:
        pxs = [a[0] for a in s_sig]
        pys = [a[1] for a in s_sig]
        pr = pearson_r(pxs, pys)
        dfc = [abs(a[0] - a[1]) for a in s_sig]
        dpc = [abs(a[2] - a[3]) for a in s_sig]
        dqc = [abs(a[4] - a[5]) for a in s_sig]
        lines.append(f"pearson_signalValue_matched\t{pr if not math.isnan(pr) else 'na'}")
        lines.append(
            f"max_abs_diff_signalValue_matched\t{max(dfc):.6g} / largest idx by diff not tracked"
        )
        lines.append(f"max_abs_diff_pValue_matched\t{max(dpc):.6g}")
        lines.append(f"max_abs_diff_qValue_matched\t{max(dqc):.6g}")
        s_sig.sort(key=lambda t: -abs(t[0] - t[1]))
        lines.append("largest5_signalValue_abs_diffs\t" + " ; ".join(f"{a[0]:.6g} vs {a[1]:.6g}" for a in s_sig[:5]))
    else:
        lines.append("pearson_signalValue_matched\tna")
        lines.append("max_abs_diff_signalValue_matched\tna")
        lines.append("max_abs_diff_pValue_matched\tna")
        lines.append("max_abs_diff_qValue_matched\tna")
        lines.append("largest5_signalValue_abs_diffs\tna")

    # rank stability: rank by pval descending per file
    rm = sorted(range(len(m)), key=lambda i: m[i].pval, reverse=True)
    rc_ = sorted(range(len(c)), key=lambda i: c[i].pval, reverse=True)
    if len(m) and len(c):
        lines.append(
            "rank_pval_top3_macs3_chroms\t" + ",".join(m[rm[i]].chrom for i in range(min(3, len(rm))))
        )
        lines.append(
            "rank_pval_top3_cpp_chroms\t" + ",".join(c[rc_[i]].chrom for i in range(min(3, len(rc_))))
        )

    # summits: compare 1bp positions
    if args.macs3_summits and args.cpp_summits:
        sm = read_smt(args.macs3_summits)
        sc = read_smt(args.cpp_summits)
        lines.append(f"summit_row_count_macs3\t{len(sm)}")
        lines.append(f"summit_row_count_cpp\t{len(sc)}")
        n_ex = 0
        n_all = min(len(sm), len(sc))
        for a, b in zip(sm[:n_all], sc[:n_all]):
            if a.chrom == b.chrom and a.start == b.start:
                n_ex += 1
        lines.append(
            f"summit_1bp_exact_same_order_rows\t{n_ex} / {n_all}" if n_all else "summit_1bp_exact_same_order_rows\t0 / 0"
        )

    out = "\n".join(lines) + "\n"
    with open(args.out_tsv, "w", encoding="utf-8") as o:
        o.write(out)
    sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
