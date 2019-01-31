#! /usr/bin/env python3
"""Prepare a phylogenetic tree for using in RASE.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import collections
import os
import sys
from ete3 import *

NW_IN_FORMAT = 1
NW_OUT_FORMAT = 1


def load_tsv_dict(tsv, col1, col2):
    """Load a dictionary from a TSV file (col1 -> col2).
    """
    d = {}
    with open(tsv) as f:
        for l in f:
            l = l.strip()
            if len(l) == 0:
                continue
            parts = l.split("\t")
            k = parts[col1]
            v = parts[col2]
            d[k] = v
    return d


def node_to_pg(node, pg_dict):
    if node.is_leaf():
        return set([pg_dict[node.name]])
    else:
        s = set()
        for n in node:
            s |= node_to_pg(n, pg_dict)
        return s


def rename_internal_nodes(tree, pg_dict):
    """Rename internal nodes (add phylogroups to the name).
    """

    numbers = collections.defaultdict(lambda: 0)

    i = 1
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        pgs = node_to_pg(node, pg_dict)
        pgs_s = "-".join(sorted(list(pgs), key=int))
        nname = "PG-{}_{}".format(pgs_s, numbers[pgs_s])
        node.name = nname
        numbers[pgs_s] += 1
    return tree


def rename_leaves(tree, rename_dict):
    """Rename leaves.
    """
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.name = rename_dict[node.name]
    return tree


def create_rase_tree(newick_in_fn, newick_out_fn, table_fn, node_col, taxid_col, pg_col):
    tree = Tree(newick_in_fn, format=NW_IN_FORMAT)

    pg_dict = load_tsv_dict(table_fn, node_col, pg_col)
    rename_dict = load_tsv_dict(table_fn, node_col, taxid_col)

    tree = rename_internal_nodes(tree, pg_dict)
    tree = rename_leaves(tree, rename_dict)

    tree.write(
        format=NW_OUT_FORMAT,
        outfile=newick_out_fn,
    )


def main():
    parser = argparse.ArgumentParser(
        description='Prepare a tree for RASE. Rename nodes and add names for internal nodes.'
    )

    parser.add_argument(
        '-i',
        '--in-newick',
        type=str,
        metavar='FILE',
        required=True,
        dest='newick_in_fn',
        help='input newick tree',
    )

    parser.add_argument(
        '-o',
        '--out-newick',
        type=str,
        metavar='FILE',
        required=True,
        dest='newick_out_fn',
        help='output newick tree',
    )

    parser.add_argument(
        '-m',
        '--metadata',
        type=str,
        metavar='FILE',
        required=True,
        dest='table_fn',
        help='TSV file',
    )

    parser.add_argument(
        '-n',
        '--node-col',
        type=int,
        metavar='INT',
        required=True,
        dest='node_col',
        help='0-based node name column id',
    )

    parser.add_argument(
        '-t',
        '--taxid-col',
        type=int,
        metavar='INT',
        required=True,
        dest='taxid_col',
        help='0-based taxid name column id (new node name)',
    )

    parser.add_argument(
        '-p',
        '--phylogroup-col',
        type=int,
        metavar='INT',
        required=True,
        dest='pg_col',
        help='0-based phylogroup column id',
    )

    args = parser.parse_args()

    create_rase_tree(
        newick_in_fn=args.newick_in_fn,
        newick_out_fn=args.newick_out_fn,
        table_fn=args.table_fn,
        node_col=args.node_col,
        taxid_col=args.taxid_col,
        pg_col=args.pg_col,
    )


if __name__ == "__main__":
    main()
