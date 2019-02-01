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

script_dir = os.path.dirname(os.path.realpath(__file__))


def debug(*vals):
    print("Debug:", *vals, file=sys.stderr)


def error(*msg, error_code=1):
    print('Rase Error:', *msg, file=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()
    sys.exit(error_code)


def rase(db, reads):

    # 0) check correctness

    # 1) decompress prophyle index

    db_dir = os.path.dirname(db)
    db_pref = os.path.basename(db)
    cmd_decompress = ['prophyle', 'decompress', db + ".tar.gz", db_dir]
    process_decompress = subprocess.Popen(cmd_decompress, shell=False)
    process_decompress.communicate()

    # 2) run prophyle classify + rase_predict.py

    cmd_classify = ['prophyle', 'classify', db, reads]
    cmd_predict = [
        'rase_predict.py',
        '-t',
        'clock',
        '-i',
        '1',
        db + "/tree.nw",
        db + '.tsv',
        '-',  #'-p', prediction/{wildcards.pref}__{wildcards.index}/
    ]

    process_classify = subprocess.Popen(cmd_classify, stdout=subprocess.PIPE, shell=False)
    process_predict = subprocess.Popen(cmd_predict, stdin=process_classify.stdout, shell=False)

    print("Running:", *cmd_classify, "|", *cmd_predict, file=sys.stderr)

    process_classify.stdout.close()
    a = process_predict.communicate()  # return ... [0]
    #print(a)


def main():
    parser = argparse.ArgumentParser(description=DESC)

    parser.add_argument('rase_db', help="prefix of the RASE database (<pref>.tar.gz, <pref>.tsv)")

    parser.add_argument('reads.fq', help="nanopore reads (- for stdin)")

    args = parser.parse_args()

    try:
        rase(db=args.rase_db, reads=args.reads)
    except KeyboardInterrupt:
        error("Keyboard interrupt")


if __name__ == "__main__":
    main()
