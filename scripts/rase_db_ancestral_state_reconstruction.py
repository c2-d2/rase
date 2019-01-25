#! /usr/bin/env python3
"""Read resistance categories and infer the missing ones (NA) using parsimony.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import collections
import copy
import csv
import ete3
import os
import re
import sys

cat_to_color = {
    "R": "red",
    "r": "violet",
    "S": "blue",
    "s": "green",
}


def dbg(*args):
    if 0:
        print(*args, file=sys.stderr)


class Parsimony:
    def __init__(self, nhx_fn, tsv_fn, antibiotics):
        self.load_tree(nhx_fn)
        self.load_res_tsv(tsv_fn)
        self.antibiotics = antibiotics

    def load_tree(self, nhx_fn):
        """Load an ete3 tree
        """
        self.tree = ete3.Tree(nhx_fn, format=1)
        for i, n in enumerate(self.tree.traverse()):
            if n.name == "":
                n.name = str("node_{}".format(i))

    def load_res_tsv(self, tsv_fn):
        self.mic_dict = collections.OrderedDict()
        self.res_dict = collections.OrderedDict()
        self.int_dict = collections.OrderedDict()

        with open(tsv_fn) as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.first_line = line.strip()
                    continue
                line = line.strip()
                taxid, mic, interval, res = line.split("\t")
                self.mic_dict[taxid] = mic
                self.int_dict[taxid] = interval
                self.res_dict[taxid] = res

    def infer_resistance_parsimony(self):
        """Infer resistance categories for isolates with NA.

        Args:
            tree: ete3 tree
            old_res_dict: taxid -> category
        """

        self._get_pars_score(self.tree)
        dbg(self.res_dict)
        dbg()
        self._propagate_categories(self.tree)
        dbg(self.res_dict)

    def plot_tree(self, pdf_fn):
        print("Plotting", pdf_fn, file=sys.stderr)
        ts = ete3.TreeStyle()
        ant = self.antibiotics.upper()
        ts.title.add_face(ete3.TextFace("{}{}".format(18 * " ", ant), fsize=60), column=0)
        ts.show_leaf_name = False

        for n in self.tree.traverse():
            taxid = n.name
            nstyle = ete3.NodeStyle()
            res_cat = self.res_dict[taxid]
            nstyle["fgcolor"] = cat_to_color[res_cat]
            nstyle["size"] = 15
            n.set_style(nstyle)

            if n.is_leaf():
                N = ete3.AttrFace("name", fsize=20)
                n.add_face(N, 0)
                #add_face_to_node(N, n, 0, position="aligned")
                #n.add_face_to_node(ete3.TextFace("XXX"))

        self.tree.render(pdf_fn, tree_style=ts)

    def _get_pars_score(self, node):
        taxid = node.name

        # The order is important here!!
        cats = ['NA', 'R', 'S']

        if node.is_leaf():
            if self.res_dict[taxid] == "R":
                score = {'R': 0, 'S': 1, 'NA': 99999}
            elif self.res_dict[taxid] == "S":
                score = {'R': 1, 'S': 0, 'NA': 99999}
            else:
                assert self.res_dict[taxid] == "NA"
                score = {'R': 0, 'S': 0, 'NA': 0}  # NA is a joker
        else:
            costs = [{'R': 0, 'S': 0, 'NA': 0} for _ in node.children]
            score = {'R': 0, 'S': 0, 'NA': 0}

            for column, n in enumerate(node.children):
                child_pars = self._get_pars_score(n)
                for to in cats:
                    candidates = []
                    for fr0m in cats:
                        if fr0m == "NA":
                            candidates.append(child_pars[fr0m])
                        elif to == "NA":
                            candidates.append(child_pars[fr0m] + 9999)
                        elif fr0m == to:
                            candidates.append(child_pars[fr0m])
                        else:
                            candidates.append(child_pars[fr0m] + 1)
                    costs[column][to] = min(candidates)

            dbg()
            dbg(taxid, [x.name for x in node.children])
            dbg("\n".join([str(c) for c in costs]))
            for c in cats:
                for j in range(len(node.children)):
                    score[c] += costs[j][c]
            dbg()

            m_score = min(score.values())
            for c in cats:
                if score[c] == m_score:
                    cat = c
                    break  # stop with first category (i.e., 'R' is preferred)

            #####
            ##### debugging - chl dataset
            #####
            #assert cat!="R"

            self.res_dict[taxid] = cat

            #assert score[2] >= max(score[:-1]), "taxid {}, score {}, children scores {}".format(taxid, score, children_scores)

        dbg("parsimony score", taxid, score)
        return score

    def _propagate_categories(self, node):
        taxid = node.name
        assert self.res_dict[taxid] != "NA"
        if node.is_leaf():
            return
        for n in node.children:
            child_taxid = n.name
            if self.res_dict[child_taxid] not in ["R", "S"]:
                self.res_dict[child_taxid] = self.res_dict[taxid].lower()
            self._propagate_categories(n)


def main(tsv_fn, nhx_fn, pdf_fn, antibiotics):
    parsimony = Parsimony(nhx_fn, tsv_fn, antibiotics)
    parsimony.infer_resistance_parsimony()
    print(parsimony.first_line)
    for taxid in parsimony.mic_dict:
        print(taxid, parsimony.mic_dict[taxid], parsimony.int_dict[taxid], parsimony.res_dict[taxid], sep="\t")
    if pdf_fn:
        parsimony.plot_tree(pdf_fn)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Read resistance categories and infer missing ones (NA) using parsimony."
    )

    parser.add_argument(
        '-i',
        '--in-tsv',
        type=str,
        metavar='str',
        required=True,
        dest='tsv_fn',
        help='Input TSV file with taxid / pseudo-mic / cat',
    )

    parser.add_argument(
        '-t',
        '--tree',
        type=str,
        metavar='str',
        required=True,
        dest='nhx_fn',
        help='Newick tree',
    )
    parser.add_argument(
        '-a',
        '--antibiotics',
        type=str,
        metavar='str',
        default='',
        dest='antibiotics',
        help='Name of the antibiotics (for titles)',
    )

    parser.add_argument(
        '-p',
        '--plot-pdf',
        type=str,
        metavar='str',
        default=None,
        dest='pdf_fn',
        help='Plot the tree with the assigned colors',
    )

    args = parser.parse_args()

    main(
        tsv_fn=args.tsv_fn,
        nhx_fn=args.nhx_fn,
        pdf_fn=args.pdf_fn,
        antibiotics=args.antibiotics,
    )
