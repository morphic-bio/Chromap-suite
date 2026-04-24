#!/usr/bin/env python3
"""
Validate paired-end mate fields in a coordinate-sorted dual (fragment-level) BAM.

For each QNAME there must be exactly two records (read1 / read2). When the mate
is mapped (no FMUNMAP), checks SAM mate chain: RNEXT/PNEXT vs the mate's
RNAME/POS (1-based), including discordant pairs on different references.

Does not replace a full synthetic alignment test; use this as a regression scan
on real Chromap dual output.
"""
from __future__ import annotations

import subprocess
import sys
from collections import defaultdict


def mate_ref(rec_rname: str, rec_rnext: str) -> str:
    if rec_rnext == "=":
        return rec_rname
    return rec_rnext


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: check_dual_bam_mate_pairs.py <file.bam>", file=sys.stderr)
        return 2
    bam = sys.argv[1]
    proc = subprocess.run(
        ["samtools", "view", bam],
        check=True,
        capture_output=True,
        text=True,
    )
    by_q: dict[str, list[list[str]]] = defaultdict(list)
    for line in proc.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 11:
            continue
        by_q[cols[0]].append(cols)

    bad = 0
    FMUNMAP = 0x8
    for qname, rows in by_q.items():
        if len(rows) != 2:
            print(
                f"WARN: expected 2 records for {qname}, got {len(rows)}",
                file=sys.stderr,
            )
            bad += 1
            continue

        def is_read1(c: list[str]) -> bool:
            return (int(c[1]) & 0x40) != 0 and (int(c[1]) & 0x80) == 0

        a, b = rows[0], rows[1]
        if is_read1(b):
            a, b = b, a
        if not is_read1(a) or is_read1(b):
            print(f"WARN: could not classify read1/read2 for {qname}", file=sys.stderr)
            bad += 1
            continue

        fa, fb = int(a[1]), int(b[1])
        # Mate unmapped on this read
        if (fa & FMUNMAP) == 0:
            exp_mate_chr = mate_ref(a[2], a[6])
            if exp_mate_chr != b[2]:
                print(
                    f"MATE_REF: {qname} read1 RNEXT/RNAME implies {exp_mate_chr}, "
                    f"mate RNAME={b[2]}",
                    file=sys.stderr,
                )
                bad += 1
            try:
                pnext = int(a[7])
                pos_b = int(b[3])
            except ValueError:
                print(f"MATE_POS: {qname} non-integer POS/PNEXT", file=sys.stderr)
                bad += 1
            else:
                if pnext != pos_b:
                    print(
                        f"MATE_POS: {qname} read1 PNEXT {pnext} != read2 POS {pos_b}",
                        file=sys.stderr,
                    )
                    bad += 1

        if (fb & FMUNMAP) == 0:
            exp_mate_chr = mate_ref(b[2], b[6])
            if exp_mate_chr != a[2]:
                print(
                    f"MATE_REF: {qname} read2 RNEXT/RNAME implies {exp_mate_chr}, "
                    f"mate RNAME={a[2]}",
                    file=sys.stderr,
                )
                bad += 1
            try:
                pnext = int(b[7])
                pos_a = int(a[3])
            except ValueError:
                print(f"MATE_POS: {qname} non-integer POS/PNEXT", file=sys.stderr)
                bad += 1
            else:
                if pnext != pos_a:
                    print(
                        f"MATE_POS: {qname} read2 PNEXT {pnext} != read1 POS {pos_a}",
                        file=sys.stderr,
                    )
                    bad += 1

    if bad:
        print(f"check_dual_bam_mate_pairs: {bad} problem(s)", file=sys.stderr)
        return 1
    print(f"check_dual_bam_mate_pairs: OK ({len(by_q)} QNAMEs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
