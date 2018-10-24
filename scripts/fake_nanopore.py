#! /usr/bin/env python3

"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import os
import sys


def message(*args):
    print(*args, file=sys.stderr)


def readfq(fp):
    ##
    ## Taken from https://github.com/lh3/readfq
    ##
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


def se_reader(reads_fn):
    with open(reads_fn) as fp:
        for name, seq, _ in readfq(fp):
            yield name, seq


def pe_reader(reads1_fn, reads2_fn):
    with open (reads1_fn) as fp1:
        with open (reads2_fn) as fp2:
            f1=readfq(fp1)
            f2=readfq(fp2)

            while 1:
                name1, seq1, _ =next(f1)
                name2, seq2, _ =next(f2)
                yield f"{name1}_{name2}", f"{seq1}n{seq2}"


def fx_stats(reads_fn, max_reads=None):
    """Calculate basic statistics for a FASTx file.

    Args:
        reads_fn(str): FASTA or FASTQ filename.
        max_reads(int): Take only first max_reads reads.

    Return
        (#reads, #bps)
    """
    if max_reads is None:
        max_reads=10**15
    reads=0
    bps=0
    with open(reads_fn) as f:
        for name, seq, qual in readfq(f):
            reads+=1
            bps+=len(seq)
            if reads >= max_reads:
                break
    return reads, bps


def first_pass(reads1_fn, reads2_fn, max_reads=None):
    """First pass through the input files.

    Args:
        reads_fn(str): 1st FASTA or FASTQ filename.
        reads_fn(str): 2nd FASTA or FASTQ filename.
        max_reads(int): Take only first x reads.

    Return
        (#reads, #bps)
    """
    reads1, bps1 = fx_stats(reads1_fn, max_reads)
    if reads2_fn:
        reads2, bps2 = fx_stats(reads1_fn, max_reads)
        assert reads1 == reads2, "The input files have different numbers of reads"
    else:
        reads2, bps2 = 0, 0
    return reads1+reads2, bps1+bps2


def fake_nanopore(reads1_fn, reads2_fn, reads_per_read, bps_per_read, number_of_reads, number_of_bps):
    """Fake nanopore reads.

    Args:
        reads1_fn(str): 1st FASTA or FASTQ filename.
        reads1_fn(str): 2nd FASTA or FASTQ filename.
        reads_per_read(int): Create reads by chaining n reads, None for no restrictions.
        bps_per_read(int): Create reads by chaining them up to x bps, None for no restrictions.
        number_of_reads(int): Take only first n reads, None for no restrictions.
        number_of_bps(): Take only first n bps, None for no restrictions.
    """


    # 1) Initialize parameters
    message("First pass through the files")
    if reads2_fn is None:
        reader=se_reader(reads1_fn)
    else:
        reader=pe_reader(reads1_fn, reads2_fn)


    total_reads, total_bps = first_pass(reads1_fn, reads2_fn, number_of_reads)

    message(f"First pass finished; {total_reads} reads and {total_bps}bps")

    if bps_per_read is None:
        bps_per_read=10**15

    if reads_per_read is None:
        reads_per_read=1

    if number_of_reads is None:
        number_of_reads=total_reads
    else:
        number_of_reads=min(total_reads, number_of_reads)

    if number_of_bps is None:
        number_of_bps=total_bps
    else:
        number_of_bps=min(reads, number_of_bps)

    # 2) Iterate over reads
    current_reads=0
    current_bps=0
    while current_reads<number_of_reads and current_bps<number_of_bps:
        names=[]
        seqs=[]
        bqs=[]

        seqs_len=0

        # 2A) Chaining
        #message("Chaining")
        while len(names)<reads_per_read \
                and seqs_len<bps_per_read:
            try:
                n,s=next(reader)
                names.append(n)
                seqs.append(s)
                seqs_len+=len(s)
            except StopIteration:
                # todo: fix end
                break

        # 2B) Updating statistics and printing
        #message("Printing")
        name=" ".join(names)
        seq="N".join(seqs)
        current_reads+=len(names)
        current_bps+=len(seq)
        fq_print_read(name, seq)


def fa_print_read(name, seq):
    print(f">{name}")
    print(seq)


def fq_print_read(name, seq, qual=None):
    print(f"@{name}")
    print(seq)
    print("+")
    if qual is None:
        print(len(seq)*"?")
    else:
        print(qual)


def ont_fake_name_gen():
    i=0
    run_id="runid=4242424242424242424242424242424242424242"
    ch="ch=42"
    while True:
        i+=1
        main_id="{:08}-4242-4242-4242-424242424242".format(i)
        read="read={}".format(i)
        start_time="start_time={}".format(fake_datetime(round(i*200/1000000)))
        yield " ".join([main_id, run_id, read, ch, start_time])


def fake_datetime(mbp):
    mbp=int(mbp)
    s=mbp%60
    mbp//=60
    m=mbp%60
    mbp//=60
    h=mbp%60
    dt="2018-01-01T{:02}:{:02}:{:02}Z".format(h, m, s)
    return dt


def main():
    desc = """\
            Create fake nanopore reads from Illumina reads or assembled contigs.
    """

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        'reads1_fn',
        metavar="reads1.fq",
        help='1st FASTA or FASTQ file',
    )

    parser.add_argument(
        'reads2_fn',
        metavar="reads2.fq",
        default=None,
        nargs='?',
        help='2nd FASTA or FASTQ file',
    )

    parser.add_argument('-r',
        dest='reads_per_read',
        metavar='INT',
        type=int,
        default=None,
        help='Create reads by chaining INT reads.',
    )

    parser.add_argument('-b',
        dest='bps_per_read',
        metavar='INT',
        type=int,
        default=None,
        help='Create reads by chaining them up to INT bps.',
    )

    parser.add_argument('-R',
        dest='number_of_reads',
        metavar='INT',
        type=int,
        default=None,
        help='Take only first INT reads',
    )

    parser.add_argument('-N',
        dest='number_of_bps',
        metavar='INT',
        type=int,
        default=None,
        help='Take only first INT bps.',
    )

    args = parser.parse_args()

    fake_nanopore(
        reads1_fn=args.reads1_fn,
        reads2_fn=args.reads2_fn,
        reads_per_read=args.reads_per_read,
        bps_per_read=args.bps_per_read,
        number_of_reads=args.number_of_reads,
        number_of_bps=args.number_of_bps,
    )


if __name__ == '__main__':
    main()
