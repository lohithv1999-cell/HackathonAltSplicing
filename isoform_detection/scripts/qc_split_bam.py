#!/usr/bin/env python3
"""Split a BAM into two halves by read name, for the split-half stability check.

The split is deterministic, so re-running gives the same two halves. Reads here are
single-end, so there are no mates to keep together, though hashing the name means
secondary records still follow their primary.
"""

import argparse
import hashlib
import pysam


def half_for(read_name):
    # md5 rather than the built-in hash(), which is randomised per process and so
    # would give a different split every run.
    return hashlib.md5(read_name.encode()).digest()[0] & 1


def main():
    parser = argparse.ArgumentParser(description="Split a BAM into two halves by read-name hash")
    parser.add_argument("-i", "--input", required=True, help="Input BAM, a pooled cell-type BAM")
    parser.add_argument("-a", "--out_a", required=True, help="Output BAM for half A")
    parser.add_argument("-b", "--out_b", required=True, help="Output BAM for half B")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="IO threads per output file [%(default)s]")
    args = parser.parse_args()

    n_a = 0
    n_b = 0

    with pysam.AlignmentFile(args.input, "rb", threads=args.threads) as infile:
        # "wb0" is uncompressed BAM. The halves get read once by the next rule and
        # then thrown away, and compressing them was taking longer than everything
        # else in the branch put together.
        with pysam.AlignmentFile(args.out_a, "wb0", template=infile, threads=args.threads) as out_a, \
             pysam.AlignmentFile(args.out_b, "wb0", template=infile, threads=args.threads) as out_b:
            write_a = out_a.write
            write_b = out_b.write
            for read in infile:
                if half_for(read.query_name) == 0:
                    write_a(read)
                    n_a += 1
                else:
                    write_b(read)
                    n_b += 1

    total = n_a + n_b
    print(f"[split] {args.input}")
    print(f"[split] half A: {n_a:,}  half B: {n_b:,}  total: {total:,}")
    if total > 0:
        print(f"[split] half A is {100 * n_a / total:.2f}% of reads (expected ~50%)")


if __name__ == "__main__":
    main()