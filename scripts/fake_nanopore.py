#! /usr/bin/env python3

"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import os
import sys

def readfq(fp):
    ## Taken from https://github.com/lh3/readfq
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def fake_nanopore_se(reads1):
    create_ont_fake_name=ont_fake_name_gen()
    for name, seq, qual in readfq(reads1):
        fake_ont_name=next(create_ont_fake_name)
        print_read(fake_ont_name, seq, qual)


def fake_nanopore_pe(reads1, reads2):
    create_ont_fake_name=ont_fake_name_gen()
    for (name1, seq1, qual1), (name2, seq2, qual2) in zip(readfq(reads1),readfq(reads1)):
        fake_ont_name=next(create_ont_fake_name)
        print_read(fake_ont_name, "{}NNN{}".format(seq1,seq2), "{}NNN{}".format(qual1,qual2))


def print_read(name, seq, qual):
        print("{}{}".format("@" if qual else ">", name))
        print(seq)
        if qual:
            print("+")
            print(qual)


def ont_fake_name_gen():
    i=0
    main_id="42424242-4242-4242-4242-424242424242"
    run_id="runid=4242424242424242424242424242424242424242"
    ch="ch=42"
    while True:
        i+=1
        read="read={}".format(i)
        start_time="start_time=2018-01-01T17:16:17Z"
        yield " ".join([main_id, run_id, read, ch, start_time])


def main():
    desc = """\
            Create fake nanopore reads from Illumina reads or assembled contigs.
    """

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        'reads1',
        metavar="reads1.fq",
        type=argparse.FileType('r'),
        help='1st FASTA or FASTQ file',
    )

    parser.add_argument(
        'reads2',
        metavar="reads2.fq",
        type=argparse.FileType('r'),
        default=None,
        nargs='?',
        help='2nd FASTA or FASTQ file',
    )

    args = parser.parse_args()
    if args.reads2 is None:
        fake_nanopore_se(args.reads1)
    else:
        fake_nanopore_pe(args.reads1, args.reads2)


if __name__ == '__main__':
    main()
