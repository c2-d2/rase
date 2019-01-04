#! /usr/bin/env python3

"""Delete a node in a Newick tree.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import os
import sys
from ete3 import *

FEATURES=[]

def delete_nodes(newick_in_fn, newick_out_fn, nodes):
    tree=Tree(newick_in_fn,format=1)

    print(nodes)
    to_delete=set(nodes)

    #for node in tree:
    #    if node.name in to_delete:
    #        del node

    for d in to_delete:
        n=tree.search_nodes(name=d)[0]
        n.delete()

    tree.write(
            format=1,
            #features=[],
            features=FEATURES,
            outfile=newick_out_fn,
            format_root_node=True,
        )


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Delete children.')

    parser.add_argument('-i', '--in-newick',
            type=str,
            metavar='FILE',
            required=True,
            dest='newick_in_fn',
            help='input newick tree',
        )

    parser.add_argument('-o', '--out-newick',
            type=str,
            metavar='FILE',
            required=True,
            dest='newick_out_fn',
            help='output newick tree',
        )

    parser.add_argument('-d', '--node-del',
            type=str,
            metavar='STR',
            nargs="+",
            dest='nodes',
            help='Nodes to delete',
            default=[],
        )


    args = parser.parse_args()

    delete_nodes(
        newick_in_fn=args.newick_in_fn,
        newick_out_fn=args.newick_out_fn,
        nodes=args.nodes,
    )
