#! /usr/bin/env python3
"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import collections
import csv
import datetime
import glob
import os
from pprint import pprint
import re
import sys

re_timestamp = re.compile(r'.*/(\d{10})\.tsv')

HEADER_PRINTED = False


class RTbl:
    """Resistance table.

    This class loads data from a file with a description of individual isolates
    and preprocesses them so that they can be accessed via its internal
    structures.

    Attributes:
        pg: taxid -> phylogroup
        pgset: phylogroup -> set of taxids
        serotype: taxid -> serotype
        st: taxid -> sequence type
        rcat: taxid -> antibiotic -> category
        measures: taxid -> measure -> value
        ants: List of antibiotics.
    """

    def __init__(self, tsv):

        self.pg = {}
        self.pgset = collections.defaultdict(set)
        self.serotype = {}
        self.st = {}

        self.rcat = collections.defaultdict(dict)

        self.measures = {}

        with open(tsv, 'r') as f:
            tsv_reader = csv.DictReader(f, delimiter='\t')

            ants = filter(lambda x: x.find("_mic") != -1, tsv_reader.fieldnames)
            self.ants = list(map(lambda x: x.replace("_mic", ""), ants))

            #print(dir(tsv_reader))
            for x in tsv_reader:
                taxid = x['taxid']

                try:
                    serotype = x['serotype']
                except KeyError:
                    serotype = "NA"

                try:
                    st = x['ST']
                except KeyError:
                    st = "NA"

                pg = x['phylogroup']
                self.pg[taxid] = pg
                self.st[taxid] = st
                self.serotype[taxid] = serotype
                self.pgset[pg].add(taxid)

                for a in self.ants:
                    cat = x[a + "_cat"].upper()
                    assert cat in set(["NA", "S", "R"])
                    self.rcat[taxid][a] = cat

        #print(self.ants)
        #print(self.pg)
        #print(self.rcat)
        #sys.exit(0)


class Tsv:
    """Assignment statistics at a time.

    This class loads prediction statistics for a given time point (relative
    similarity to individual samples).

    Attributes:
        rtbl: Resistance table
        phylogroups: Sorted list of phylogroups
        filename: TSV filename
        timestamp: Unix timestamp of that sequencing time
        datetime: The corresponding datetime
        measures: Measures  (taxid -> measure -> value)
    """

    #############
    # Interface #
    #############

    def __init__(self, tsv, rtbl):
        self.rtbl = rtbl
        self.phylogroups = sorted(self.rtbl.pg.keys())
        self.filename = tsv

        m = re_timestamp.match(tsv)
        if m:
            self.timestamp = int(m.group(1))
            dt = datetime.datetime.fromtimestamp(self.timestamp)
            self.datetime = dt.strftime('%Y-%m-%d %H:%M:%S')
        else:
            self.timestamp = 'NA'
            self.datetime = 'NA'

        self.load_measures(tsv)

    def load_measures(self, tsv):
        """Load measures from the provided TSV.

        Set self.measures (taxid -> measure -> value).
        """

        print("Loading", tsv, file=sys.stderr)

        with open(tsv, 'r') as f:
            tsv_reader = csv.DictReader(f, delimiter='\t')
            rows = []
            for r in tsv_reader:
                #if r["taxid"]=="_unassigned_":
                #    continue

                rd = {}
                ###
                rd['taxid'] = r['taxid']
                # count
                rd['count'] = float(r['count'])
                # ln
                rd['ln'] = float(r['ln'])
                # h1
                rd['h1'] = float(r['h1'])
                # c1
                rd['c1'] = float(r['c1'])
                rows.append(rd)

        rows.sort(key=lambda x: x['h1'], reverse=True)
        self.measures = {r['taxid']: r for r in rows}

    def predict(self, measure):
        """Predict.
        """

        summary = collections.OrderedDict()

        ppgs = self._predict_pg(measure=measure)

        sorted_pgs = list(ppgs)
        predicted_pg = sorted_pgs[0]
        pg1_taxid, pg1_measmax = ppgs[sorted_pgs[0]]
        pg2_taxid, pg2_measmax = ppgs[sorted_pgs[1]]

        predicted_taxid = pg1_taxid
        predicted_serotype = self.rtbl.serotype[predicted_taxid]
        predicted_st = self.rtbl.st[predicted_taxid]

        ####

        summary['datetime'] = self.datetime
        summary['read count'] = int(self.cumul_count())
        summary['read len'] = int(self.cumul_ln())
        summary['used bases'] = int(self.cumul_c1())

        summary['PG1'] = sorted_pgs[0]
        summary['PG1_meas'] = round(pg1_measmax)
        #summary['PG1meas']='%s' % float('%.3g' % pg1_measmax)
        summary['PG2'] = sorted_pgs[1]
        summary['PG2_meas'] = round(pg2_measmax)

        summary['taxid'] = predicted_taxid
        summary['serotype'] = predicted_serotype
        summary['ST'] = predicted_st

        if pg1_measmax == 0:
            summary['PG1'] = "NA"
            summary['taxid'] = "NA"
            summary['serotype'] = "NA"
            summary['ST'] = "NA"
        if pg2_measmax == 0:
            summary['PG2'] = "NA"
        if pg1_measmax > 0:
            summary['PG_score'] = 2 * (round(pg1_measmax / (pg1_measmax + pg2_measmax) * 1000) / 1000) - 1
        else:
            summary['PG_score'] = 0

        for ant in self.rtbl.ants:
            pres = self._predict_resistance(predicted_pg, ant)

            # ant category
            cat_col = ant.upper() + "_cat"
            if pg1_measmax > 0:
                predict_cat = self.rtbl.rcat[predicted_taxid][ant]
            else:
                predict_cat = "NA"
            summary[cat_col] = predict_cat

            # todo: add a comment that we take the category of the best isolate; not the same as in the plots

            # susc score
            score_col = ant.upper() + "_susc_score"
            try:
                s_meas = pres['S'][1]
                r_meas = pres['R'][1]
                if r_meas + s_meas > 0:
                    susc_score = round(1000 * s_meas / (r_meas + s_meas)) / 1000
                else:
                    susc_score = 0
            except KeyError:
                # computing susc score fails
                if predict_cat == 'R':
                    susc_score = 0.0
                elif predict_cat == 'S':
                    susc_score = 1.0
                elif predict_cat == 'NA' and pg1_measmax == 0:
                    susc_score = 0.0
                elif predict_cat == 'NA':
                    susc_score = 'NA'

            summary[score_col] = susc_score

        self.summary = summary

    def print(self):
        """Print.
        """

        global HEADER_PRINTED

        summary = self.summary

        keys = list(summary)
        values = [summary[k] for k in summary]
        if not HEADER_PRINTED:
            print(*keys, sep="\t")
            HEADER_PRINTED = True
        print(*values, sep="\t")

    #############
    # Auxiliary #
    #############

    def _predict_pg(self, measure="h1"):
        """Predict phylogroup.

        Returns:
            SortedDict: pg -> (taxid, val)
        """

        d = collections.defaultdict(lambda: [None, -1])

        for taxid in self.measures:
            if taxid == "_unassigned_":
                continue
            pg = self.rtbl.pg[taxid]
            val = self.measures[taxid][measure]

            if val > d[pg][1]:
                d[pg] = taxid, val

        l = list(d.items())
        l.sort(key=lambda x: x[1][1], reverse=True)
        return collections.OrderedDict(l)

    def _predict_resistance(self, pg, ant, measure="h1"):
        """Predict resistance for a given phylogroup.

        ...for given antibiotics and with respect to given measure.

        Returns:
            category -> (taxid, val)
        """

        d = collections.defaultdict(lambda: (None, -1))

        for taxid in self.rtbl.pgset[pg]:
            val = self.measures[taxid][measure]
            cat = self.rtbl.rcat[taxid][ant]

            if val > d[cat][1]:
                d[cat] = taxid, val

        return dict(d)

    def cumul_count(self):
        return sum([self.measures[taxid]['count'] for taxid in self.measures])

    def cumul_ln(self):
        return sum([self.measures[taxid]['ln'] for taxid in self.measures])

    def cumul_c1(self):
        return sum([self.measures[taxid]['c1'] for taxid in self.measures])


def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        'res_tbl',
        type=str,
        metavar='res_tbl.tsv',
    )

    parser.add_argument(
        'file',
        type=str,
        metavar='file.tsv',
        nargs='+',
    )

    args = parser.parse_args()

    rtbl = RTbl(args.res_tbl)

    #summary(args.file)
    for f in args.file:
        tsv = Tsv(f, rtbl)
        tsv.predict("h1")
        tsv.print()


if __name__ == "__main__":
    main()
