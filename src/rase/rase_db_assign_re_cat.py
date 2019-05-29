#! /usr/bin/env python3
"""Read MICs from a TSV file, parse the intervals and assign resistance categories (R / S / NA).

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import os
import sys
import re
import csv

re_number = re.compile(r'^([0-9]*\.{0,1}[0-9]*)$')
re_greater = re.compile(r'^(?:≥|>|>=)([0-9]*\.{0,1}[0-9]*)$')
re_lesser = re.compile(r'^(?:≤|<|<=)([0-9]*\.{0,1}[0-9]*)$')


def pseudo_mic_to_interval(pseudo_mic):
    """Parse MIC string and return the corresponding interval.
    """

    pseudo_mic = pseudo_mic.strip()
    pseudo_mic = pseudo_mic.replace(" ", "")
    pseudo_mic = pseudo_mic.split("/")[0]

    if pseudo_mic in ["-", "", "NA"]:
        return (0, float("+inf"))

    m = re_number.match(pseudo_mic)
    if m:
        number = float(m.group(1))
        return (number, number)

    m = re_greater.match(pseudo_mic)
    if m:
        number = float(m.group(1))
        return (number, float("+inf"))

    m = re_lesser.match(pseudo_mic)
    if m:
        number = float(m.group(1))
        return (0, number)

    print(
        f"Warning: MIC string '{pseudo_mic}' could not be parsed. If this represents a value, a new regular expression should be added.",
        file=sys.stderr)
    return (0, float("+inf"))


def interval_to_cat(interval, breakpoint):
    """Compare an MIC interval with a breakpoint and assign a resistance category.
    """
    x, y = interval
    assert x <= y

    if x >= breakpoint:
        return "R"
    elif y < breakpoint:
        return "S"
    else:
        return "NA"


def assign_cat(
        tsv_in_fn,
        ant,
        breakpoint,
        miccol,
        taxidcol,
):

    with open(tsv_in_fn) as tsv_fo:
        #reader = csv.DictReader(tsv_fo, delimiter='\t')

        print("\t".join(["taxid", ant + "_mic", ant + "_int", ant + "_cat"]))

        for i, row in enumerate(tsv_fo):
            if i == 0:
                continue

            parts = row.split("\t")
            #print(row)
            taxid = parts[taxidcol].strip()
            pseudo_mic = parts[miccol].strip()
            interval = pseudo_mic_to_interval(pseudo_mic)
            cat = interval_to_cat(interval, breakpoint)

            if pseudo_mic == "":
                pseudo_mic = "NA"

            line = "\t".join(
                [taxid, pseudo_mic, "{}-{}".format(*interval), cat])
            line = line.replace("\tinf", "\tInf")  # correction for R
            line = line.replace("-inf", "-Inf")  # correction for R
            print(line)


def main():
    parser = argparse.ArgumentParser(
        description=
        "Read MICs from a TSV file and assign resistance categories (R / S / NA)."
    )

    parser.add_argument(
        '-i',
        '--in-tsv',
        type=str,
        metavar='str',
        required=True,
        dest='tsv_in_fn',
        help='Input TSV file with taxid / pseudo-mic / cat',
    )

    parser.add_argument(
        '-t',
        '--taxid-col',
        type=int,
        metavar='int',
        required=True,
        dest='taxidcol',
        help='0-based taxid column ID',
    )

    parser.add_argument(
        '-m',
        '--mic-col',
        type=int,
        metavar='str',
        required=True,
        dest='miccol',
        help='0-based MIC column ID',
    )

    parser.add_argument(
        '-b',
        '--breakpoint',
        type=float,
        metavar='float',
        required=True,
        dest='breakpoint',
        help='Breakpoint for antibiotic resistance',
    )

    parser.add_argument(
        '-a',
        '--antibiotic',
        type=str,
        metavar='str',
        default='ant',
        dest='ant',
        help='Antibiotic acronym (to print) [ant]',
    )

    args = parser.parse_args()

    assign_cat(
        tsv_in_fn=args.tsv_in_fn,
        ant=args.ant,
        breakpoint=args.breakpoint,
        miccol=args.miccol,
        taxidcol=args.taxidcol,
    )


if __name__ == "__main__":
    main()
