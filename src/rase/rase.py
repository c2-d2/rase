#! /usr/bin/env python3
"""rase

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import os
import subprocess
import sys

sys.path.append(os.path.dirname(__file__))
import version

PROGRAM = 'rase'
VERSION = version.VERSION
DESC = 'RASE - rapid prediction of antibiotic resistance using lineage calling'

GITDIR = os.path.basename(sys.argv[0])[-3:] == ".py"
if GITDIR:
    C_D = os.path.abspath(os.path.dirname(sys.argv[0]))
else:
    C_D = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

RASE = os.path.join(C_D, "rase")
RASE_PREDICT = os.path.join(C_D, "rase_predict.py")


def debug(*vals):
    print("Debug:", *vals, file=sys.stderr)


def error(*msg, error_code=1):
    print('Rase Error:', *msg, file=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()
    sys.exit(error_code)


def rase(db, reads_fn, bam_fn):

    # 0) check correctness

    db_tsv = db + ".tsv"
    db_tgz = db + ".tar.gz"
    for x in (db_tsv, db_tgz):
        if not os.path.isfile(db_tsv):
            error("File '{}' could not be found".format(x))
    if reads_fn != "-" and not os.path.isfile(reads_fn):
        error("File '{}' could not be found".format(reads_fn))

    # 1) decompress prophyle index

    db_dir = os.path.dirname(db)
    db_pref = os.path.basename(db)
    cmd_decompress = ['prophyle', 'decompress', db_tgz, db_dir]
    process_decompress = subprocess.Popen(cmd_decompress, shell=False)
    process_decompress.communicate()

    # 2) run prophyle classify + rase_predict.py

    cmd_classify = ['prophyle', 'classify', db, reads_fn]
    cmd_predict = [
        RASE_PREDICT,
        '-t',
        'clock',
        '-s',
        '1',
        db + "/tree.nw",
        db + '.tsv',
        '-',  #'-p', prediction/{wildcards.pref}__{wildcards.index}/
    ]
    if bam_fn is not None:
        cmd_predict.append(bam_fn[0])

    # todo: maybe use subprocess.check_output ?
    process_classify = subprocess.Popen(cmd_classify, stdout=subprocess.PIPE, shell=False)
    process_predict = subprocess.Popen(cmd_predict, stdin=process_classify.stdout, shell=False)

    print("Running:", *cmd_classify, "|", *cmd_predict, file=sys.stderr)

    process_classify.stdout.close()
    a = process_predict.communicate()  # return ... [0]
    #print(a)


def main():
    parser = argparse.ArgumentParser(description=DESC)

    parser.add_argument('rase_db', metavar='db', help="prefix of the RASE database (<pref>.tar.gz, <pref>.tsv)")

    parser.add_argument('reads', metavar='reads.fq', help="nanopore reads (- for stdin)")

    parser.add_argument(
        'out_bam', metavar='assignments.bam', help="computed ProPhyle RASE assignments", nargs="*", default=[]
    )

    args = parser.parse_args()

    if args.out_bam:
        out_bam = args.out_bam
    else:
        out_bam = None

    try:
        rase(db=args.rase_db, reads_fn=args.reads, bam_fn=out_bam)
    except KeyboardInterrupt:
        error("Keyboard interrupt")


if __name__ == "__main__":
    main()
