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
DESC = ''

script_dir = os.path.dirname(os.path.realpath(__file__))


def run(*args):
    print("Running: '{}'".format(' '.join(args)), file=sys.stderr)
    subprocess.call(args)


def rase(db, reads):
    db_dir = os.path.dirname(db)
    db_pref = os.path.basename(db)
    run(*['cd', db_dir, '&&', 'prophyle', 'decompress', db_pref+".tar.gz"])
    run(*['prophyle', 'classify', db, reads,
          '|', script_dir+'/../scripts/rase_predict.py', db+"/tree.nw", db+'.tsv', '-', #'-p', prediction/{wildcards.pref}__{wildcards.index}/
          ])
          #> {params.tsv1}

def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        'rase_db'
    )

    parser.add_argument(
        'reads',
    )

    args = parser.parse_args()
    rase(db=args.rase_db, reads=args.reads)


if __name__ == "__main__":
    main()
