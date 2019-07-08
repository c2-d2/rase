#! /usr/bin/env python3
"""Compute statistics about resistance from the main DB table.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import collections
import csv
import re

re_cat = re.compile("(.*)_cat")
cats = ['S', 'R', 's', 'r', 'Rr', 'Ss']


def compute_stats(fn, s):
    if s:
        cats = ['S', 'R', 's', 'r', 'Rr', 'Ss']
    else:
        cats = ['Rr', 'Ss']

    stats = collections.defaultdict(lambda: collections.Counter())
    ants = {}
    lineages = set()
    lineage_count = collections.Counter()

    with open(fn) as f:
        tsv_reader = csv.DictReader(f, delimiter='\t')
        for r in tsv_reader:
            try:
                lineage = r["lineage"]
            except KeyError:
                lineage = r["pg"]
            except KeyError:
                lineage = r["phylogroup"]
            lineage_count[lineage] += 1
            lineages.add(lineage)

            for k in r:
                m = re_cat.match(k)
                if m:
                    ant = m.group(1)
                    cat = r[k]
                    ants[ant] = ''

                    stats[lineage][(ant, cat)] += 1

    ants_sorted = sorted(ants.keys(), key=lambda s: s.lower())
    cats_sorted = sorted(cats, key=lambda s: s.lower())

    parts = ['lineage', 'count'] + [
        '{}_{}'.format(ant, cat) for ant in ants_sorted for cat in cats_sorted
    ]
    print(*parts, sep="\t")

    for lineage in sorted(lineages, key=lambda x: int(x)):
        parts = [lineage, lineage_count[lineage]]
        for ant in ants_sorted:
            for catgr in cats_sorted:
                parts.append(
                    sum([stats[lineage][(ant, cat)] for cat in catgr]))

        print(*parts, sep="\t")

    parts = [
        "sum",
        sum([lineage_count[lineage] for lineage in lineages]),
    ]
    for ant in ants_sorted:
        for catgr in cats_sorted:
            parts.append(
                sum([
                    sum([stats[lineage][(ant, cat)] for cat in catgr])
                    for lineage in lineages
                ]))

    print(*parts, sep="\t")


def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        '-s',
        action="store_false",
        dest="s",
        help='Shorter version of the outputs (only _Rr and _Ss categories)',
    )

    parser.add_argument(
        'fn',
        type=str,
        metavar='res_cat.tsv',
        help='File with resistance categories',
    )

    args = parser.parse_args()

    compute_stats(args.fn, args.s)


if __name__ == "__main__":
    main()
