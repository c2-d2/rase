#! /usr/bin/env python3
"""Compute statistics about resistance from the main DB table.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import collections
import csv
import itertools
import os
import re
import sys

re_cat = re.compile("(.*)_cat")
cats = ['S', 'R', 's', 'r', 'Rr', 'Ss']


def compute_stats(fn):
    stats = collections.defaultdict(lambda: collections.Counter())
    ants = {}
    pgs = set()
    pg_count = collections.Counter()

    with open(fn) as f:
        tsv_reader = csv.DictReader(f, delimiter='\t')
        for r in tsv_reader:
            pg = r["phylogroup"]
            pg_count[pg] += 1
            pgs.add(pg)

            for k in r:
                m = re_cat.match(k)
                if m:
                    ant = m.group(1)
                    cat = r[k]
                    ants[ant] = ''

                    stats[pg][(ant, cat)] += 1

            #print(pg,r)

    ants_sorted = sorted(ants.keys(), key=lambda s: s.lower())
    cats_sorted = sorted(cats, key=lambda s: s.lower())

    parts = ['pg', 'count'] + ['{}_{}'.format(ant, cat) for ant in ants_sorted for cat in cats_sorted]
    print(*parts, sep="\t")

    for pg in sorted(pgs, key=lambda x: int(x)):
        parts = [pg, pg_count[pg]]
        for ant in ants_sorted:
            for catgr in cats_sorted:
                parts.append(sum([stats[pg][(ant, cat)] for cat in catgr]))

        print(*parts, sep="\t")


def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        'fn',
        type=str,
        metavar='res_cat.tsv',
        help='File with resistance categories',
    )

    args = parser.parse_args()

    compute_stats(args.fn)


if __name__ == "__main__":
    main()
