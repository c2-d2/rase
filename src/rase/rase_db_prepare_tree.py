#! /usr/bin/env python3
"""Prepare phylogenetic tree for use in RASE.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import collections
import hashlib
import os
import pprint as ppprint
import sys
from ete3 import *

NW_IN_FORMAT = 1
NW_OUT_FORMAT = 1


def test_file(fn):
    try:
        with open(fn) as _:
            pass
    except FileNotFoundError:
        print(f"File '{fn}' doesn't exists")
        sys.exit(1)


def pprint(x, max_len=100):
    s = ppprint.pformat(x)
    if len(s) > 100:
        print(s[:100] + "...\n")
    else:
        print(s)


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


def node_to_lineage(node, lineage_dict):
    if node.is_leaf():
        try:
            return set([lineage_dict[node.name]])
        except KeyError:
            print(
                f"Renaming error: Unknown isolate '{node.name}'",
                file=sys.stderr)
            sys.exit(42)
    else:
        s = set()
        for n in node:
            s |= node_to_lineage(n, lineage_dict)
        return s


def _sorting_key(s):
    """Numerical/lexicographical sorting.
    """
    try:
        val = int(s)
        norm = "{:010}".format(val)
    except ValueError:
        norm = s
    return norm


def rename_internal_nodes(tree, lineage_dict):
    """Rename internal nodes (add phylogroups to the name).
    """

    numbers = collections.defaultdict(lambda: 0)

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        lineages = node_to_lineage(node, lineage_dict)
        lineages_s = "-".join(sorted(list(lineages), key=_sorting_key))
        # if the name is too long, make a bibtex-like etc form and attach hash
        if len(lineages_s)<40:
            lineages_ss = lineages_s
        else:
            # short "fake" name
            x = "-".join(sorted(list(lineages), key=_sorting_key)[:3])
            h = hashlib.md5(lineages_ss.encode('utf-8')).hexdigest()[:10]
            lineages_ss = f"{x}-etc--{h}"
        nname = "L{}_{}".format(lineages_ss, numbers[lineages_s])
        node.name = nname
        numbers[lineages_s] += 1
    return tree


def rename_leaves(tree, rename_dict):
    """Rename leaves.
    """
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.name = rename_dict[node.name]
    return tree


def prepare_rase_tree(newick_in_fn, newick_out_fn, table_fn, node_col,
                      taxid_col, lineage_col):
    print("\n1) Testing files.\n")
    test_file(newick_in_fn)
    test_file(table_fn)
    tree = Tree(newick_in_fn, format=NW_IN_FORMAT)

    print("\n2) Building phylogroup dictionary:\n")
    lineage_dict = load_tsv_dict(table_fn, node_col, lineage_col)
    pprint(lineage_dict)

    print("\n3) Building renaming dictionary:\n")
    rename_dict = load_tsv_dict(table_fn, node_col, taxid_col)
    pprint(rename_dict)

    print("\n4) Renaming nodes.\n")
    tree = rename_internal_nodes(tree, lineage_dict)
    tree = rename_leaves(tree, rename_dict)

    tree.write(
        format=NW_OUT_FORMAT,
        outfile=newick_out_fn,
    )


def main():
    parser = argparse.ArgumentParser(
        description=
        'Prepare a tree for RASE. Rename nodes and add names for internal nodes.'
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
        dest='lineage_col',
        help='0-based phylogroup column id',
    )

    args = parser.parse_args()

    prepare_rase_tree(
        newick_in_fn=args.newick_in_fn,
        newick_out_fn=args.newick_out_fn,
        table_fn=args.table_fn,
        node_col=args.node_col,
        taxid_col=args.taxid_col,
        lineage_col=args.lineage_col,
    )


if __name__ == "__main__":
    main()
