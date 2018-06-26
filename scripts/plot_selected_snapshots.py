#! /usr/bin/env python3

import argparse
import collections
import datetime
import glob
import os
import re
import subprocess
import sys

re_fn=re.compile(r'.*?(\d+)\.tsv')


def extract_ts(fn):
    """Extract timestamp from a file name
    """
    m=re_fn.match(fn)
    ts=m.group(1)
    return ts


def secs_to_text(cs):
    """Convert time in seconds to a human readable string.
    """
    v=collections.OrderedDict()
    v['s']=cs%60
    cs//=60
    v['m']=cs%60
    cs//=60
    v['h']=cs%24
    cs//=24
    v['d']=cs
    parts=[]
    for k in v:
        if v[k]!=0:
            parts.append("{}{}".format(v[k],k))
    return "_".join(parts[::-1])


def plot_snapshots(res_table, directory, indexes, outprefix, time0):
    plotting_script=os.path.join(os.path.dirname(os.path.realpath(__file__)), "plot_snapshot.R")
    plots=get_plot_info(directory, indexes, time0)
    for fn, s, t in plots:
        cmd=[plotting_script, res_table, fn, "{}{}.pdf".format(outprefix, t)]
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)


def get_plot_info(directory, indexes, time0):
    fns=glob.glob(os.path.join(directory, "*.tsv"))
    tsvs_abs={
            fn: int(extract_ts(fn))
            for fn in fns
            }
    min_ts=min(tsvs_abs.values())
    tsvs_rel=[
            [k, v-min_ts+time0]
            for (k,v) in tsvs_abs.items()
            ]
    tsvs_rel.sort(key=lambda x: x[1])

    tr=[]
    for i in indexes:
        fn, sec = tsvs_rel[i]
        tr.append((fn, sec, secs_to_text(sec)))

    return tr


def main():
    parser = argparse.ArgumentParser(description="Plot selected snapshots and use human readable names (e.g., 5m, 1h, etc.)")

    parser.add_argument('-f',
            type=int,
            default=60,
            metavar='int',
            dest='time0',
            help='First timestamp [60]',
        )

    parser.add_argument('res_table',
            type=str,
            metavar='res_table.tsv',
            help='Resistance table',
        )

    parser.add_argument('dir',
            type=str,
            metavar='dir',
            help='Snapshot directory',
        )

    parser.add_argument('indexes',
            type=int,
            metavar='index',
            nargs='+',
            help='0-based index of a snapshot (e.g., 0, -2, 6)',
        )

    parser.add_argument('outpref',
            type=str,
            metavar='outpref',
            help='Output prefix',
        )

    args = parser.parse_args()

    plot_snapshots(args.res_table, args.dir, args.indexes, args.outpref, args.time0)


if __name__ == "__main__":
    main()

