#! /usr/bin/env python3

"""rase

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import os
import sys

sys.path.append(os.path.dirname(__file__))
import version

PROGRAM = 'rase'
VERSION = version.VERSION
DESC = ''

def rase(db, reads):
    cmd='prophyle decompress "{}"'.format(db)
    cmd="""
    prophyle decompress
    prophyle classify {index} {reads}
    ../scripts/rase_predict.py {params.tree} {params.metadata} {params.bam} -p prediction/{wildcards.pref}__{wildcards.index}/ > {params.tsv1}
    """


def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        'rase_db'
    )

    parser.add_argument(
        'reads',
    )

    #parser.add_argument(
    #    '-v',
    #    '--version',
    #    action='version',
    #    version='{} {}'.format(PROGRAM, VERSION),
    #)

    args = parser.parse_args()
    rase(db=args.rase_db, reads=args.reads)


if __name__ == "__main__":
    main()
