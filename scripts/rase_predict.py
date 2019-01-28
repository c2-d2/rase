#! /usr/bin/env python3
"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

# ./scripts/rase_predict.py ~/github/my/rase-predict/database/spneumoniae_sparc.k18/tree.nw ~/github/my/rase-predict/prediction/sp10_norwich_P33.filtered__spneumoniae_sparc.k18.bam ~/github/my/rase-predict/database/spneumoniae_sparc.k18.tsv | tL

import argparse
import collections
import csv
import datetime
import ete3
import itertools
import glob
import json
import os
import pysam
import re
import sys
import warnings

FAKE_ISOLATE_UNASSIGNED = "_unassigned_"
HEADER_PRINTED = False
re_timestamp = re.compile(r'.*/(\d{10})\.tsv')


def timestamp_from_qname(qname):
    return int(qname.partition("_")[0])


def format_time(seconds):
    minutes = seconds // 60
    hours = minutes // 60
    minutes -= 60 * hours
    return "{}h{}m".format(hours, minutes)


class Runner:
    def __init__(self, metadata_fn, tree_fn, bam_fn, pref, delta, first_read_delay):
        self.metadata = RaseMetadataTable(metadata_fn)
        self.stats = Stats(tree_fn)
        self.predict = Predict(self.metadata)
        self.rase_bam_reader = RaseBamReader(bam_fn)
        self.pref = pref
        self.delta = delta
        self.first_read_delay = first_read_delay

    def run(self):
        # 1) set the initial window:
        #     [t0, t0+delta), where t0=time_of_first_read-first_read_delay
        t0 = self.rase_bam_reader.t1 - self.first_read_delay
        current_window = [t0, t0 + self.delta]  # [x, y)

        if self.pref is not None:
            f = open("{}/{}.tsv".format(self.pref, current_window[1]), mode="w")

        # 2) iterate through individual reads, and update and print statistics
        for read_stats in self.rase_bam_reader:
            assert len(read_stats) > 0
            read_timestamp = timestamp_from_qname(read_stats[0]["qname"])

            # do we have to shift the window?
            if read_timestamp >= current_window[1]:
                if self.pref is not None:
                    self.stats.print(file=f)
                self.predict.predict(self.stats)
                self.predict.print()
                if self.pref is not None:
                    f.close()
                while read_timestamp >= current_window[1]:
                    current_window[0] += self.delta
                    current_window[1] += self.delta
                if self.pref is not None:
                    f = open("{}/{}.tsv".format(self.pref, current_window[1]), mode="w")

                last_time = format_time(current_window[0] - t0)
                print(
                    "Time t={}: {} reads and {} non-propagated ({} propagated) assignments processed.".format(
                        last_time, self.stats.nb_assigned_reads, self.stats.nb_nonprop_asgs, self.stats.nb_asgs
                    ), file=sys.stderr
                )

            self.stats.update_from_read(read_stats)

        if self.pref is not None:
            self.stats.print(file=f)
            f.close()

        # todo: check if the last version is printed, if not, uncomment
        #with open(args.tsv, mode="w") as f:
        #    stats.print(file=f)


class Predict:
    """Predicting from assignment statistics at a time.

    This class loads prediction statistics for a given time point (relative
    similarity to individual samples).

    Attributes:
        rtbl: Resistance table
        phylogroups: Sorted list of phylogroups
        filename: TSV filename
        timestamp: Unix timestamp of that sequencing time
        datetime: The corresponding datetime
        weight: taxid -> weight
    """

    def __init__(self, rtbl):
        self.rtbl = rtbl
        self.phylogroups = sorted(self.rtbl.pgset.keys())
        self.summary = collections.OrderedDict()
        self.taxids = sorted(self.rtbl.pg.keys())

    def predict(self, stats):
        """Predict.
        """

        ppgs = self._predict_pg(stats.stats_h1)
        sorted_pgs = list(ppgs)
        predicted_pg = sorted_pgs[0]
        pg1_taxid, pg1_measmax = ppgs[sorted_pgs[0]]
        pg2_taxid, pg2_measmax = ppgs[sorted_pgs[1]]
        predicted_taxid = pg1_taxid

        #predicted_serotype = self.rtbl.serotype[predicted_taxid]
        #predicted_st = self.rtbl.st[predicted_taxid]

        ####

        current_dt, _, _ = str(datetime.datetime.now()).partition(".")


        self.summary['datetime'] = current_dt #"now" #todo: self.datetime
        #self.summary['read count'] = int(self.cumul_count())
        #self.summary['read len'] = int(self.cumul_ln())
        self.summary['PG1'] = sorted_pgs[0]
        self.summary['PG1_meas'] = round(pg1_measmax)
        self.summary['PG2'] = sorted_pgs[1]
        self.summary['PG2_meas'] = round(pg2_measmax)
        self.summary['taxid'] = predicted_taxid
        #self.summary['serotype'] = predicted_serotype
        #self.summary['ST'] = predicted_st

        if pg1_measmax == 0:
            self.summary['PG1'] = "NA"
            self.summary['taxid'] = "NA"
            #self.summary['serotype'] = "NA"
            #self.summary['ST'] = "NA"
        if pg2_measmax == 0:
            self.summary['PG2'] = "NA"
            # todo: PG2 taxid
        if pg1_measmax > 0:
            self.summary['PG_score'] = 2 * (round(pg1_measmax / (pg1_measmax + pg2_measmax) * 1000) / 1000) - 1
        else:
            self.summary['PG_score'] = 0

#        for ant in self.rtbl.ants:
#            pres = self._predict_resistance(predicted_pg, ant)
#
#            # ant category
#            cat_col = ant.upper() + "_cat"
#            if pg1_measmax > 0:
#                predict_cat = self.rtbl.rcat[predicted_taxid][ant]
#            else:
#                predict_cat = "NA"
#            self.summary[cat_col] = predict_cat
#
#            # todo: add a comment that we take the category of the best isolate; not the same as in the plots
#
#            # susc score
#            score_col = ant.upper() + "_susc"
#            try:
#                s_meas = pres['S'][1]
#                r_meas = pres['R'][1]
#                if r_meas + s_meas > 0:
#                    susc_score = round(1000 * s_meas / (r_meas + s_meas)) / 1000
#                else:
#                    susc_score = 0
#            except KeyError:
#                # computing susc score fails
#                if predict_cat == 'R':
#                    susc_score = 0.0
#                elif predict_cat == 'S':
#                    susc_score = 1.0
#                elif predict_cat == 'NA' and pg1_measmax == 0:
#                    susc_score = 0.0
#                elif predict_cat == 'NA':
#                    susc_score = 'NA'
#
#            self.summary[score_col] = susc_score

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

    def _predict_pg(self, weight):
        """Predict phylogroup.

        Returns:
            SortedDict: pg -> (taxid, val)
        """

        d = collections.defaultdict(lambda: [None, -1])

        for taxid in self.taxids:
            if taxid == "_unassigned_":
                continue
            pg = self.rtbl.pg[taxid]
            val = weight[taxid]

            if val > d[pg][1]:
                d[pg] = taxid, val

        l = list(d.items())
        l.sort(key=lambda x: x[1][1], reverse=True)
        return collections.OrderedDict(l)

    def _predict_resistance(self, pg, ant):
        """Predict resistance for a given phylogroup.

        ...for given antibiotics and with respect to h1.

        Returns:
            category -> (taxid, val)
        """

        d = collections.defaultdict(lambda: (None, -1))

        for taxid in self.rtbl.pgset[pg]:
            val = self.weight[taxid]
            cat = self.rtbl.rcat[taxid][ant]

            if val > d[cat][1]:
                d[cat] = taxid, val

        return dict(d)

    #def cumul_count(self):
    #    return sum([self.measures[taxid]['count'] for taxid in self.measures])
    #
    #def cumul_ln(self):
    #    return sum([self.measures[taxid]['ln'] for taxid in self.measures])


class Stats:
    """Statistics for RASE predictions.

    The attributes contain all necessary information that are necessary for
    predicting in RASE.

    Params:
        tree_fn (str): Filename of the Newick tree.

    Attributes:
        nb_reads (int): Number of processed reads.
        nb_nonprop_asgs (int): Number of processed alignments (before propagation).
        nb_asgs (int): Number of processed alignments (after propagation).

        cumul_qlen (float): Cumulative weighted read length.
        cumul_h1 (float): Cumulative weighted hit count.

        stats_ct (dict): isolate -> number of processed reads
        stats_h1 (dict): isolate -> squared h1
        stats_ln (dict): isolate -> weighted qlen
    """

    def __init__(self, tree_fn):
        self._tree = ete3.Tree(tree_fn, format=1)
        self._isolates = sorted([isolate.name for isolate in self._tree])
        self._descending_isolates_d = self._precompute_descendants(self._tree)

        # stats for assigned reads
        self.nb_assigned_reads = 0
        self.nb_nonprop_asgs = 0
        self.nb_asgs = 0

        # stats for unassigned reads
        self.nb_unassigned_reads = 0

        # cumulative statistics for assigned reads
        self.cumul_h1 = 0.0
        self.cumul_ln = 0.0

        # statistics for individual isolates, "_unassigned_" for unassigned
        self.stats_ct = collections.defaultdict(lambda: 0.0)
        self.stats_h1 = collections.defaultdict(lambda: 0.0)
        self.stats_ln = collections.defaultdict(lambda: 0.0)

    @staticmethod
    def _precompute_descendants(tree):
        descending_leaves = {}
        for root in list(tree.traverse()) + [tree]:
            descending_leaves[root.name] = set([isolate.name for isolate in root])
        return descending_leaves

    def _descending_isolates(self, *nnames):
        isolates = set()
        for nname in nnames:
            isolates |= self._descending_isolates_d[nname]
        return sorted(isolates)

    def update_from_read(self, asgs):
        """Update statistics from assignments of a single read.

        Params:
            asgs (dict): Assignments.
        """

        assert len(asgs) > 0, "No assignments provided"

        asg0 = asgs[0]
        is_assigned = asg0["assigned"]

        # is read is assigned to at least 1 node?
        if is_assigned:
            h1 = asg0['h1']
            ln = asg0['ln']

            # 1) update cumulative stats
            self.nb_assigned_reads += 1
            self.nb_nonprop_asgs += len(asgs)

            self.cumul_h1 += h1
            self.cumul_ln += ln

            # 2) update stats for individual isolates
            nnames = [asg["rname"] for asg in asgs]
            bestmatching_isolates = self._descending_isolates(*nnames)
            l = len(bestmatching_isolates)
            delta_ct = 1.0 / l
            delta_h1 = h1 / l
            delta_ln = ln / l
            for isolate in bestmatching_isolates:
                self.stats_ct[isolate] += delta_ct
                self.stats_h1[isolate] += delta_h1
                self.stats_ln[isolate] += delta_ln
        else:
            if len(asgs) != 1:
                warnings.warn("A single read shouldn't be reported as unassigned multiple times ({})".format(asgs))
            self.nb_unassigned_reads += 1
            self.update_isolate_stats([FAKE_ISOLATE_UNASSIGNED], h1=0, ln=asg0["ln"], l=1)

    def print(self, file):
        """Print statistics to a file.

        Args:
            file (file): Output file.
        """
        print("taxid", "count", "count_norm", "ln", "ln_norm", "h1", "h1_norm", sep="\t", file=file)
        table = []
        for isolate in self._isolates + [FAKE_ISOLATE_UNASSIGNED]:
            table.append(
                [
                    isolate,
                    self.stats_ct[isolate],
                    1.0 * self.stats_ct[isolate] / self.nb_assigned_reads if self.nb_assigned_reads != 0 else 0,
                    self.stats_ln[isolate],
                    self.stats_ln[isolate] / self.cumul_ln if self.cumul_ln != 0 else 0,
                    self.stats_h1[isolate],
                    self.stats_h1[isolate] / self.cumul_h1 if self.cumul_h1 != 0 else 0,
                ]
            )

        table.sort(key=lambda x: x[5], reverse=True)

        for x in table:
            print(*x, sep="\t", file=file)


class RaseMetadataTable:
    """RASE metadata table.

    This class loads data from a file with a description of individual isolates
    and preprocesses them so that they can be accessed via its internal
    structures.

    Attributes:
        pg: taxid -> phylogroup
        pgset: phylogroup -> set of taxids
        serotype: taxid -> serotype
        st: taxid -> sequence type
        rcat: taxid -> antibiotic -> category
        weight: taxid -> weight
        ants: List of antibiotics.
    """

    def __init__(self, tsv):

        self.pg = {}
        self.pgset = collections.defaultdict(set)
        self.serotype = {}
        self.st = {}

        self.rcat = collections.defaultdict(dict)

        self.weight = {}

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


class RaseBamReader:
    """Iterator over all assignments of individual reads in a BAM/RASE file.

    Assumes a non-empty BAM file.

    Params:
        bam_fn (str): BAM file name.

    Attributes:
        t1 (int): Timestamp of first read
    """

    def __init__(self, bam_fn):
        self.assignment_reader = SingleAssignmentReader(bam_fn)
        self._buffer = []
        self._finished = False
        self._load_assignment()
        self._extract_t1()

    def __iter__(self):
        return self

    def __next__(self):
        """Retrieve set of assignments of the next read.
        """

        if self._finished:
            raise StopIteration

        while len(self._buffer) < 2 or self._buffer[-1]["qname"] == self._buffer[-2]["qname"]:
            try:
                self._load_assignment()
            except StopIteration:
                self._finished = True
                buffer = self._buffer
                self._buffer = []
                return buffer

        buffer = self._buffer[:-1]
        self._buffer = self._buffer[-1:]
        return buffer

    def _load_assignment(self):
        asg = next(self.assignment_reader)
        self._buffer.append(asg)

    def _extract_t1(self):
        self.t1 = timestamp_from_qname(self._buffer[0]["qname"])


class SingleAssignmentReader:
    """Iterator over individual assignments in BAM/RASE.

    Assumes that it is possible to infer read lengths (either
    from the ln tag, base sequence or cigars).

    Params:
        bam_fn (str): BAM file name.

    Attributes:
        samfile (pysam.AlignmentFile): PySAM file object.
    """

    def __init__(self, bam_fn):
        self.samfile = pysam.AlignmentFile(bam_fn, "rb")
        self.alignments_iter = self.samfile.fetch(until_eof=True)

        self.qname = None
        self.last_qname = None
        self.read_ln = None

    def __iter__(self):
        return self

    def __next__(self):
        """Get next assignment.

        Returns:
            assignment (dict): A dict with the following keys: "rname", "qname", "qlen, "assigned", "h1".
        """

        # 1) acquire a new alignment from SAM
        try:
            alignment = next(self.alignments_iter)
        except StopIteration:
            raise StopIteration

        # 2) re-compute basic read-related variables if necessary
        self.qname = alignment.qname
        if self.qname != self.last_qname:
            try:
                self.read_ln = alignment.get_tag("ln")
            except KeyError:  # a ln tag is not present
                if alignment.seq != "*":
                    self.read_ln = len(alignment.seq)
                else:
                    self.read_ln = read.infer_read_length()
        self.last_qname = self.qname

        # 3) infer assignment-related variables
        if not alignment.is_unmapped:
            asg = {
                "rname": alignment.reference_name,
                "qname": alignment.qname,
                "assigned": True,
                "ln": self.read_ln,
                "h1": alignment.get_tag("h1"),
            }
        else:
            asg = {
                "rname": None,
                "qname": alignment.qname,
                "assigned": False,
                "ln": self.read_ln,
                "h1": None,
            }

        return asg


def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        'tree',
        type=str,
        metavar='<tree.nw>',
    )

    parser.add_argument(
        'bam',
        type=argparse.FileType('r'),
        metavar='<assignments.bam>',
    )

    parser.add_argument(
        'metadata',
        type=str,
        metavar='metadata.tsv',
    )

    parser.add_argument('-p', type=str, dest='pref', metavar='STR', help="Output dir for samplings", default=None)

    parser.add_argument(
        '-i',
        type=int,
        dest='delta',
        metavar='INT',
        help='sampling interval (in seconds) [300]',
        default=300,
    )

    parser.add_argument(
        '-f',
        type=int,
        dest='first_read_delay',
        metavar='INT',
        help='delay of the first read [60]',
        default=60,
    )

    args = parser.parse_args()

    r = Runner(
        metadata_fn=args.metadata, tree_fn=args.tree, bam_fn=args.bam, pref=args.pref, delta=args.delta,
        first_read_delay=args.first_read_delay
    )
    r.run()


if __name__ == "__main__":
    main()
