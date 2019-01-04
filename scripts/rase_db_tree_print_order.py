#! /usr/bin/env python3
""" Print order of nodes

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import os
import sys
from ete3 import *


def print_order(newick_in_fn):
    tree = Tree(newick_in_fn, format=1)

    leaves = []

    for i, leaf in enumerate(tree.iter_leaf_names(), 1):
        leaves.append((leaf, i))

    leaves.sort(key=lambda x: x[0])

    print("taxid\torder")
    for l, i in leaves:
        print("{}\t{}".format(l, i))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Print order of leaves')

    parser.add_argument(
        '-i',
        '--in-newick',
        type=str,
        metavar='FILE',
        required=True,
        dest='newick_in_fn',
        help='input newick tree',
    )

    args = parser.parse_args()

    print_order(newick_in_fn=args.newick_in_fn, )
