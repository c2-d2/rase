#! /usr/bin/env python3
"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import datetime
import functools
import gzip
import json
import os
import re
import sys

#rre=re.compile(r'^>([^\s]+).*start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})')
#rre=re.compile(r'^>([0-9a-f\-]+) runid=([0-9a-f]) read=(\d+) ch=(\d+) start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})')

#rre=re.compile(r'''^>([0-9a-f\-]+) runid=([0-9a-f]+) read=(\d+) ch=(\d+) start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})''')#, flags=re.X)
#rre=re.compile(r'''^>([^\s]+)(.*) \s
#    start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})''', flags=re.X)

re_minion = re.compile(
    r'''^
        runid=([0-9a-f]+) \s
        read=(\d+) \s
        ch=(\d+) \s
        start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})
    ''', flags=re.X
)


def parse_header_line(hl):
    """Take a FA/FQ header line and parse all Minion fields.

    Args:
        hl (str): Header line.

    Returns:
        (uuid, runid, read, ch, date, time)
    """
    m = re_minion.match(hl)
    gg = m.groups()
    assert gg is not None, "Header line doesn't match the regular expression ('{}').".format(hl)
    return gg


def timestamp_from_datetime(date, time):
    """Get a timestamp from a date and time.

    Args:
        date (str): Date ('YY-MM-DD').
        time (str): Time ('hh:mm:ss').

    Returns:
        timestamp (int): Unix timestamp.
    """
    date_time = date + " " + time
    b = datetime.datetime.strptime(date_time, "%Y-%m-%d %H:%M:%S")
    ts = int(b.timestamp())
    return ts


def readfq(fp):
    """
    Parse FASTA/FASTQ.

    Original version: https://github.com/lh3/readfq/blob/master/readfq.py
    Changes: reads also comments, returns (name, comment, seq, quals)
    """
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        name, _, comment = last[1:].partition(" ")
        seqs = []
        last = None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, comment, ''.join(seqs), None  # yield a fasta record
            if not last: break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, comment, seq, ''.join(seqs)
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, comment, seq, None  # yield a fasta record instead
                break


def transf(fq):
    d = {}

    if fq[-3:] == ".gz":
        open_fn = functools.partial(gzip.open, mode="rt")
    else:
        open_fn = open

    with open_fn(fq) as f:
        for name, comment, seq, qual in readfq(f):
            #print(comment, file=sys.stderr)
            uuid = name
            runid, read, ch, date, time = parse_header_line(comment)
            ts = timestamp_from_datetime(date, time)
            print(">" + "{}_EX{}_RD{}_CH{}".format(ts, runid[:9], uuid.partition("-")[0], ch))
            print(seq)

        return

        for i, x in enumerate(f):
            if i % 4 == 0:
                uuid, runid, read, ch, date, time = parse_header_line(x)
                ts = timestamp_from_datetime(date, time)
                print("@" + "{}_EX{}_RD{}_CH{}".format(ts, runid[:9], uuid.partition("-")[0], ch))
            else:
                print(x, end="")


def main():
    parser = argparse.ArgumentParser(description="Reformat Minion read names in a FQ file.")

    parser.add_argument(
        'fq',
        type=str,
        metavar='reads.fq',
    )

    args = parser.parse_args()

    transf(args.fq)


if __name__ == "__main__":
    main()
