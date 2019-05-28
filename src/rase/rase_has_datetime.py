#! /usr/bin/env python3
"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT

Check if first line of a given file has a MinION header. If yes, return 0 and 1 otherwise.
"""

import argparse
import functools
import gzip
import re
import sys

# original
re_minion = re.compile(
    r'''
        runid=([0-9a-f]+) \s
        read=(\d+) \s
        ch=(\d+) \s
        start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})
    ''', flags=re.X
)

# postprocessed by RASE
re_minion_postprocessed = re.compile(
    r'''
        \d+_EX[0-9a-f]+_RD[0-9a-f]+_CH\d+
    ''',flags=re.X
    )   


def has_datetime(line):
    """Check if line contains a MinION header.

    Args:
        line (str): Header line.

    Returns:
        (uuid, runid, read, ch, date, time)
    """
    m = re_minion.search(line)
    if m:
        return True
    m = re_minion_postprocessed.search(line)
    if m:
        return True
    else:
        return False


def check_first_line(fn):
    print(f"Checking datetime in '{fn}'", file=sys.stderr)
    if fn[-3:] == ".gz":
        open_fn = functools.partial(gzip.open, mode="rt")
    else:
        open_fn = open

    with open_fn(fn) as f:
        for line in f:
            print("   checking first line: '{}...'".format(line[:30]), file=sys.stderr)
            return has_datetime(line)
        return False


def main():
    parser = argparse.ArgumentParser(description="Check whether a given FASTA/FASTQ file has MinION datetimes.")

    parser.add_argument(
        'fq',
        type=str,
        metavar='reads.fq',
    )

    args = parser.parse_args()
    r = check_first_line(args.fq)
    if r:
        print("   datetime found", file=sys.stderr)
        sys.exit(0)
    else:
        print("   datetime not found", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
