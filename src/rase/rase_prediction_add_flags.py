#! /usr/bin/env python3
"""
Add flags to RASE predictions.

S: Final state stabilized.
D: Final state detected.
L: Final state lost.

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import collections
import csv
import re


def load_tsv_dl(tsv_fn):
    """Load TSV as a list of dictionaries (1 dict per 1 line).
    """
    with open(tsv_fn) as f:
        tr = [x for x in csv.DictReader(f, delimiter="\t")]
    return tr


def get_cols(tsv_dl):
    """Get columns from a list of dictionaries.
    """
    d = tsv_dl[0]
    return [x for x in d.keys()]


def extract_flag_cols(cols):
    """Extract columns for which flags are to be computed.
    """
    res_flagging = [
        re.compile(r"^bm$"),
        re.compile(r"^serotype$"),
        re.compile(r"^bm_\S+$"),
        re.compile(r"^pg$"),
        re.compile(r"^.*_pred$"),
    ]
    flag_cols = []
    for k in cols:
        for re_flagging in res_flagging:
            if re_flagging.match(k):
                flag_cols.append(k)
                break

    return flag_cols


def get_stabilization_point(tsv_dl, key):
    """Find points of stabilization.
    """
    tsv_dl_ = [collections.defaultdict(lambda: "")
               ] + tsv_dl  # probably to support empty files?
    final_value = tsv_dl_[-1][key]
    for i in range(len(tsv_dl_) - 1, -1, -1):
        if tsv_dl_[i][key] != final_value:
            break
    return i


def add_flags_line():
    for k in flag_cols:
        if i == stabilization_points[k]:
            flags.append('S:{}'.format(k))
            assert rec[k] != prev_rec[k]
        if rec[k] != prev_rec[k]:
            if rec[k] == last_rec[k]:
                # detected
                flags.append('D:{}'.format(k))
            elif prev_rec[k] == last_rec[k]:
                # lost successfuly detected
                flags.append('L:{}'.format(k))
    flags.sort()
    flags_str = str(flags).replace("'", "")
    print(*rec.values(), flags_str, sep="\t")


def add_flags(tsv_fn):
    """Read RASE prediction output and add flags.
    """
    tsv_dl = load_tsv_dl(tsv_fn)
    cols = get_cols(tsv_dl)
    flag_cols = extract_flag_cols(cols)
    last_rec = tsv_dl[-1]
    stabilization_points = {
        col: get_stabilization_point(tsv_dl, col)
        for col in flag_cols
    }

    # add flags
    prev_rec = collections.defaultdict(lambda: "")
    for i, rec in enumerate(tsv_dl):
        if i == 0:
            print(*rec.keys(), "flags", sep="\t")
            flag_cols = extract_flag_cols(rec.keys())
        flags = []
        for k in flag_cols:
            if i == stabilization_points[k]:
                flags.append('S:{}'.format(k))
                assert rec[k] != prev_rec[k]
            if rec[k] != prev_rec[k]:
                if rec[k] == last_rec[k]:
                    # detected
                    flags.append('D:{}'.format(k))
                elif prev_rec[k] == last_rec[k]:
                    # lost successfuly detected
                    flags.append('L:{}'.format(k))
        flags.sort()
        flags_str = str(flags).replace("'", "")
        print(*rec.values(), flags_str, sep="\t")
        prev_rec = rec


def main():
    parser = argparse.ArgumentParser(
        description=
        "Add flags to predictions (R=beginning of a run, S=stabilized).")

    parser.add_argument(
        'file',
        type=str,
        metavar='predictions.tsv',
        help='',
    )

    args = parser.parse_args()
    add_flags(args.file)


if __name__ == "__main__":
    main()
