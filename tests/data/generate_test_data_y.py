#!/usr/bin/env python3
"""Generate deterministic test data with chr1 + chrY for BAM/Y split tests."""

import argparse
import os
import random


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(base, "N") for base in reversed(seq))


def rand_seq(length, rng):
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_fastq(path, reads, read_len):
    qual = "I" * read_len
    with open(path, "w") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate paired-end test data with chr1 + chrY"
    )
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    rng = random.Random(1337)
    read_len = 100
    contig_len = 400

    chr1 = rand_seq(contig_len, rng)
    chry = rand_seq(contig_len, rng)

    ref_path = os.path.join(args.output, "test_ref_y.fa")
    with open(ref_path, "w") as ref_fh:
        ref_fh.write(">chr1\n")
        ref_fh.write(chr1 + "\n")
        ref_fh.write(">chrY\n")
        ref_fh.write(chry + "\n")

    pairs = [
        ("chr1_read1", chr1, 0, 160),
        ("chr1_read2", chr1, 40, 200),
        ("chrY_read1", chry, 0, 160),
        ("chrY_read2", chry, 40, 200),
    ]

    r1_reads = []
    r2_reads = []
    for name, contig, r1_start, r2_start in pairs:
        r1_seq = contig[r1_start : r1_start + read_len]
        r2_seq = reverse_complement(contig[r2_start : r2_start + read_len])
        r1_reads.append((f"{name}/1", r1_seq))
        r2_reads.append((f"{name}/2", r2_seq))

    r1_path = os.path.join(args.output, "test_pe_y_R1.fq")
    r2_path = os.path.join(args.output, "test_pe_y_R2.fq")
    write_fastq(r1_path, r1_reads, read_len)
    write_fastq(r2_path, r2_reads, read_len)

    print(f"Generated test data in {args.output}")
    print(f"  - {ref_path}")
    print(f"  - {r1_path}")
    print(f"  - {r2_path}")


if __name__ == "__main__":
    main()
