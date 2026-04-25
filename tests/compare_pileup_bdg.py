#!/usr/bin/env python3
"""Compare two bedGraph files (chrom, start, end, value). Diagnostic parity only.

Computes total signal, covered bp (value > eps), alignment on merged constant
segments, and **genomic-span-weighted** Pearson/Spearman/MAE (each segment
weighted by its bp length). Chromosomes present in only one track contribute
segments with value 0 on the missing side. High-signal BED overlap is unchanged.
"""
from __future__ import annotations

import argparse
import bisect
import math
import sys
from typing import List, Tuple


Interval = Tuple[str, int, int, float]  # chrom, start, end, value (half-open)


def load_bdg(path: str, eps: float) -> List[Interval]:
    out: List[Interval] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, s, e, v = parts[0], parts[1], parts[2], parts[3]
            try:
                a, b = int(s), int(e)
                val = float(v)
            except ValueError:
                continue
            if b <= a:
                continue
            if abs(val) < eps:
                continue
            out.append((chrom, a, b, val))
    out.sort(key=lambda x: (x[0], x[1], x[2]))
    return out


def load_bdg_all_rows(path: str) -> List[Interval]:
    """All bedGraph intervals including value ~0 (for exact file parity)."""
    out: List[Interval] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, s, e, v = parts[0], parts[1], parts[2], parts[3]
            try:
                a, b = int(s), int(e)
                val = float(v)
            except ValueError:
                continue
            if b <= a:
                continue
            out.append((chrom, a, b, val))
    out.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    return out


def bdg_rows_identical(a: List[Interval], b: List[Interval], tol: float) -> bool:
    if len(a) != len(b):
        return False
    for (c1, s1, e1, v1), (c2, s2, e2, v2) in zip(a, b):
        if c1 != c2 or s1 != s2 or e1 != e2:
            return False
        if abs(v1 - v2) > tol:
            return False
    return True


def total_signal_and_covered(ivs: List[Interval], eps: float) -> Tuple[float, int]:
    sig = 0.0
    cov = 0
    for _c, a, b, v in ivs:
        if abs(v) < eps:
            continue
        w = b - a
        sig += w * v
        cov += w
    return sig, cov


def merge_breakpoints(
    chrom: str, one: List[Tuple[int, int, float]], two: List[Tuple[int, int, float]]
) -> List[Tuple[int, int]]:
    pts = set()
    for a, b, _ in one + two:
        pts.add(a)
        pts.add(b)
    return sorted(pts)


def segments_for_chrom(ivs: List[Interval], chrom: str) -> List[Tuple[int, int, float]]:
    return [(a, b, v) for c, a, b, v in ivs if c == chrom]


def prep_nonoverlapping(
    segs: List[Tuple[int, int, float]],
) -> Tuple[List[Tuple[int, int, float]], List[int]]:
    segs = sorted(segs, key=lambda t: (t[0], t[1]))
    starts = [s for s, _e, _v in segs]
    return segs, starts


def value_at_nonoverlap(pos: int, segs: List[Tuple[int, int, float]], starts: List[int]) -> float:
    i = bisect.bisect_right(starts, pos) - 1
    if i < 0:
        return 0.0
    s, e, v = segs[i]
    if s <= pos < e:
        return v
    return 0.0


def align_segments(
    a_iv: List[Interval], b_iv: List[Interval]
) -> Tuple[
    List[float], List[float], List[float], List[str], List[int], List[int]
]:
    """Piecewise-constant (x,y) pairs on merged breakpoints; weights = segment bp.

    Also returns per-segment (chrom, start, end) for diagnostic intervals.
    """
    chroms = sorted(set(c for c, _, _, _ in a_iv) | set(c for c, _, _, _ in b_iv))
    xs: List[float] = []
    ys: List[float] = []
    wts: List[float] = []
    segs_ch: List[str] = []
    segs_s: List[int] = []
    segs_e: List[int] = []
    for ch in chroms:
        sa = segments_for_chrom(a_iv, ch)
        sb = segments_for_chrom(b_iv, ch)
        if not sa and not sb:
            continue
        if not sa:
            sb_o, _sb_s = prep_nonoverlapping(sb)
            for a, b, v in sb_o:
                span = float(b - a)
                if span <= 0:
                    continue
                xs.append(0.0)
                ys.append(v)
                wts.append(span)
                segs_ch.append(ch)
                segs_s.append(int(a))
                segs_e.append(int(b))
            continue
        if not sb:
            sa_o, _sa_s = prep_nonoverlapping(sa)
            for a, b, v in sa_o:
                span = float(b - a)
                if span <= 0:
                    continue
                xs.append(v)
                ys.append(0.0)
                wts.append(span)
                segs_ch.append(ch)
                segs_s.append(int(a))
                segs_e.append(int(b))
            continue
        sa_o, sa_s = prep_nonoverlapping(sa)
        sb_o, sb_s = prep_nonoverlapping(sb)
        pts = merge_breakpoints(ch, sa, sb)
        for i in range(len(pts) - 1):
            p0, p1 = pts[i], pts[i + 1]
            if p1 <= p0:
                continue
            mid = (p0 + p1 - 1) // 2
            if mid < p0:
                mid = p0
            va = value_at_nonoverlap(mid, sa_o, sa_s)
            vb = value_at_nonoverlap(mid, sb_o, sb_s)
            w = float(p1 - p0)
            if w <= 0:
                continue
            xs.append(va)
            ys.append(vb)
            wts.append(w)
            segs_ch.append(ch)
            segs_s.append(int(p0))
            segs_e.append(int(p1))
    return xs, ys, wts, segs_ch, segs_s, segs_e


def wpearson(x: List[float], y: List[float], w: List[float]) -> float:
    if not x:
        return float("nan")
    sw = sum(w)
    if sw <= 0:
        return float("nan")
    mx = sum(ww * xx for xx, ww in zip(x, w)) / sw
    my = sum(ww * yy for yy, ww in zip(y, w)) / sw
    vx = sum(ww * (xx - mx) ** 2 for xx, ww in zip(x, w))
    vy = sum(ww * (yy - my) ** 2 for yy, ww in zip(y, w))
    cxy = sum(ww * (xx - mx) * (yy - my) for xx, yy, ww in zip(x, y, w))
    if vx <= 0 or vy <= 0:
        return float("nan")
    return cxy / math.sqrt(vx * vy)


def rankdata(values: List[float]) -> List[float]:
    n = len(values)
    idx = list(range(n))
    idx.sort(key=lambda i: values[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[idx[j + 1]] == values[idx[i]]:
            j += 1
        avg = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[idx[k]] = avg
        i = j + 1
    return ranks


def wspearman(x: List[float], y: List[float], w: List[float]) -> float:
    if not x:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    return wpearson(rx, ry, w)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("a_bdg", help="bedGraph A (e.g. MACS3 pileup)")
    ap.add_argument("b_bdg", help="bedGraph B (e.g. chromap_callpeaks --pileup-bdg)")
    ap.add_argument("--eps", type=float, default=1e-9)
    ap.add_argument(
        "--exact-row-tol",
        type=float,
        default=1e-4,
        help="max |Δvalue| for bdg_rows_identical_sorted (half-open coords must match)",
    )
    ap.add_argument("--high-frac", type=float, default=0.01, help="top fraction of signal mass for high-signal beds")
    ap.add_argument(
        "--high-prefix",
        default="",
        help="If set, write high-signal BEDs as <prefix>_a.bed and <prefix>_b.bed",
    )
    ap.add_argument(
        "--largest-n-diff",
        type=int,
        default=0,
        help="if >0, print top N merged segments by bp*abs(a-b) and optional --largest-out TSV",
    )
    ap.add_argument(
        "--largest-out",
        default="",
        help="Write largest-diff intervals (tab) when --largest-n-diff>0",
    )
    args = ap.parse_args()

    a_iv = load_bdg(args.a_bdg, args.eps)
    b_iv = load_bdg(args.b_bdg, args.eps)
    sa, ca = total_signal_and_covered(a_iv, args.eps)
    sb, cb = total_signal_and_covered(b_iv, args.eps)

    xs, ys, wts, seg_ch, seg_s, seg_e = align_segments(a_iv, b_iv)
    sw = sum(wts) if wts else 0.0
    pear = wpearson(xs, ys, wts) if xs else float("nan")
    spear = wspearman(xs, ys, wts) if xs else float("nan")
    if xs and sw > 0:
        mae = (
            sum(w * abs(xs[i] - ys[i]) for i, w in enumerate(wts)) / sw
        )
        mre = (
            sum(
                w
                * (
                    abs(xs[i] - ys[i])
                    / max(1e-12, max(abs(xs[i]), abs(ys[i])))
                )
                for i, w in enumerate(wts)
            )
            / sw
        )
    else:
        mae = float("nan")
        mre = float("nan")
    max_abs = max((abs(xs[i] - ys[i]) for i in range(len(xs))), default=float("nan"))

    # Exclude bp where both piecewise values are ~0 (large inter-interval gaps
    # dominate sw and deflate MAE / change correlation).
    nz_idx = [
        i
        for i in range(len(xs))
        if abs(xs[i]) >= args.eps or abs(ys[i]) >= args.eps
    ]
    mx = [xs[i] for i in nz_idx]
    my = [ys[i] for i in nz_idx]
    mw = [wts[i] for i in nz_idx]
    sw_nz = sum(mw) if mw else 0.0
    pear_nz = wpearson(mx, my, mw) if mx else float("nan")
    spear_nz = wspearman(mx, my, mw) if mx else float("nan")
    if mx and sw_nz > 0:
        mae_nz = sum(mw[i] * abs(mx[i] - my[i]) for i in range(len(mx))) / sw_nz
        mre_nz = (
            sum(
                mw[i]
                * (
                    abs(mx[i] - my[i])
                    / max(1e-12, max(abs(mx[i]), abs(my[i])))
                )
                for i in range(len(mx))
            )
            / sw_nz
        )
    else:
        mae_nz = float("nan")
        mre_nz = float("nan")
    max_abs_nz = max(
        (abs(mx[i] - my[i]) for i in range(len(mx))), default=float("nan")
    )

    print("metric\tvalue")
    print(f"file_a\t{args.a_bdg}")
    print(f"file_b\t{args.b_bdg}")
    print(f"total_signal_a\t{sa:.12g}")
    print(f"total_signal_b\t{sb:.12g}")
    print(f"total_signal_abs_diff\t{abs(sa - sb):.12g}")
    print(f"total_signal_rel_diff\t{(abs(sa - sb) / max(1e-12, max(abs(sa), abs(sb)))):.12g}")
    print(f"covered_bp_a\t{ca}")
    print(f"covered_bp_b\t{cb}")
    print(f"segment_count_merged\t{len(xs)}")
    print(f"genomic_span_weighted_bp\t{sw:.12g}")
    print(f"pearson_bp_weighted\t{pear:.12g}")
    print(f"spearman_bp_weighted\t{spear:.12g}")
    print(f"mae_bp_weighted\t{mae:.12g}")
    print(f"max_abs_diff_segment\t{max_abs:.12g}")
    print(f"mean_rel_err_bp_weighted\t{mre:.12g}")
    print(
        "note_bp_weighted\t"
        "Includes inter-interval gaps where both tracks are 0; "
        "see *_either_nonempty rows for signal-bearing bp only."
    )
    print(f"segment_count_either_nonempty\t{len(nz_idx)}")
    print(f"genomic_span_either_nonempty_bp\t{sw_nz:.12g}")
    print(f"pearson_bp_weighted_either_nonempty\t{pear_nz:.12g}")
    print(f"spearman_bp_weighted_either_nonempty\t{spear_nz:.12g}")
    print(f"mae_bp_weighted_either_nonempty\t{mae_nz:.12g}")
    print(f"max_abs_diff_segment_either_nonempty\t{max_abs_nz:.12g}")
    print(f"mean_rel_err_bp_weighted_either_nonempty\t{mre_nz:.12g}")

    ra = load_bdg_all_rows(args.a_bdg)
    rb = load_bdg_all_rows(args.b_bdg)
    ex = bdg_rows_identical(ra, rb, args.exact_row_tol)
    print(f"bdg_row_count_a\t{len(ra)}")
    print(f"bdg_row_count_b\t{len(rb)}")
    print(f"bdg_rows_identical_sorted\t{ex}")

    if args.largest_n_diff > 0 and xs and seg_ch:
        err_mass = [wts[i] * abs(xs[i] - ys[i]) for i in range(len(xs))]
        order = sorted(range(len(err_mass)), key=lambda i: err_mass[i], reverse=True)
        take = order[: min(args.largest_n_diff, len(order))]
        rows_out: List[str] = []
        for rnk, j in enumerate(take, start=1):
            print(
                f"largest_diff_{rnk}\t{seg_ch[j]}\t{seg_s[j]}\t{seg_e[j]}\t"
                f"{xs[j]:.12g}\t{ys[j]:.12g}\t{err_mass[j]:.12g}\t{wts[j]:.12g}"
            )
            rows_out.append(
                f"{seg_ch[j]}\t{seg_s[j]}\t{seg_e[j]}\t{xs[j]:.12g}\t{ys[j]:.12g}\t{err_mass[j]:.12g}\t{wts[j]:.12g}\n"
            )
        if args.largest_out:
            with open(args.largest_out, "w", encoding="utf-8") as fo:
                fo.write("chrom\tstart\tend\tval_a\tval_b\tabsdiff_mass_bp\twidth_bp\n")
                for line in rows_out:
                    fo.write(line)
            print(f"largest_diff_tsv\t{args.largest_out}")

    # High-signal intervals: greedy take intervals by mass until frac reached
    def high_bed(ivs: List[Interval], frac: float, suffix: str) -> str:
        items = []
        for c, a, b, v in ivs:
            w = (b - a) * abs(v)
            if w <= 0:
                continue
            items.append((w, c, a, b))
        items.sort(reverse=True)
        target = frac * sum(w for w, _, _, _ in items)
        acc = 0.0
        rows: List[Tuple[str, int, int]] = []
        for w, c, a, b in items:
            rows.append((c, a, b))
            acc += w
            if acc >= target:
                break
        path = (
            f"{args.high_prefix}_{suffix}.bed"
            if args.high_prefix
            else f"/tmp/compare_pileup_high_{suffix}.bed"
        )
        rows.sort(key=lambda x: (x[0], x[1], x[2]))
        with open(path, "w", encoding="utf-8") as out:
            for c, a, b in rows:
                out.write(f"{c}\t{a}\t{b}\n")
        print(f"high_signal_bed_{suffix}\t{path}")
        return path

    if a_iv and b_iv and args.high_frac > 0:
        high_bed(a_iv, args.high_frac, "a")
        high_bed(b_iv, args.high_frac, "b")
        print("high_signal_note\tbedtools jaccard on high_signal_bed paths (diagnostic)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
