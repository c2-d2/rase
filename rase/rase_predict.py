#! /usr/bin/env python3
"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

# ./scripts/rase_predict.py ~/github/my/rase-predict/database/spneumoniae_sparc.k18/tree.nw ~/github/my/rase-predict/database/spneumoniae_sparc.k18.tsv  ~/github/my/rase-predict/prediction/sp10_norwich_P33.filtered__spneumoniae_sparc.k18.bam | tL

# todo: add non-susc threshold as a param

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

import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

FAKE_ISOLATE_UNASSIGNED = "_unassigned_"
HEADER_PRINTED = False
re_timestamp = re.compile(r'.*/(\d{10})\.tsv')
re_timestamped_read = re.compile(r'\d+_.*')


def debug(*vals):
    print("Debug:", *vals, file=sys.stderr)


def error(*msg, error_code=1):
    print('Rase(predict) Error:', *msg, file=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()
    sys.exit(error_code)


def timestamp_from_qname(qname):
    if re_timestamped_read.match(qname):
        return int(qname.partition("_")[0])
    else:
        return None


def current_datetime():
    current_dt, _, _ = str(datetime.datetime.now()).partition(".")
    return current_dt


def current_timestamp():
    dto = datetime.datetime.now()
    return round(dto.timestamp())


def format_time(seconds):
    minutes = seconds // 60
    hours = minutes // 60
    minutes -= 60 * hours
    return "{}h{}m".format(hours, minutes)


def format_floats(*values, digits=3):
    form_values = []
    fstr = "{0:." + str(digits) + "f}"
    for x in values:
        if isinstance(x, float):
            form_values.append(fstr.format(round(x, 3)))
        else:
            form_values.append(x)
    return form_values


class Runner:
    def __init__(self, metadata_fn, tree_fn, bam_fn, out_bam_fn, pref, mode, delta, first_read_delay):
        self.mode = mode
        self.metadata = RaseMetadataTable(metadata_fn)
        self.stats = Stats(tree_fn, self.metadata)
        self.predict = Predict(self.metadata)
        self.rase_bam_reader = RaseBamReader(bam_fn, out_bam_fn)
        self.pref = pref
        self.delta = delta
        self.first_read_delay = first_read_delay

    def run(self):
        # 1) set the initial window:
        #     [t0, t0+delta), where t0=time_of_first_read-first_read_delay
        if self.mode == "clock":
            t0 = current_timestamp() - self.first_read_delay
        elif self.mode == "name":
            t0 = self.rase_bam_reader.t1 - self.first_read_delay
        else:
            assert (1 == 2)
        current_window = [t0, t0 + self.delta]  # [x, y)

        if self.pref is not None:
            f = open("{}/{}.tsv".format(self.pref, current_window[1]), mode="w")

        # 2) iterate through individual reads, and update and print statistics
        for read_stats in self.rase_bam_reader:
            assert len(read_stats) > 0
            if self.mode == "clock":
                read_timestamp = current_timestamp()
            elif self.mode == "name":
                read_timestamp = timestamp_from_qname(read_stats[0]["qname"])
            else:
                assert (1 == 2)

            # do we have to shift the window?
            if read_timestamp >= current_window[1]:
                if self.pref is not None:
                    self.stats.print(file=f)
                self.predict.predict(self.stats)
                self.predict.print(read_timestamp)
                if self.pref is not None:
                    f.close()
                while read_timestamp >= current_window[1]:
                    current_window[0] += self.delta
                    current_window[1] += self.delta
                if self.pref is not None:
                    f = open("{}/{}.tsv".format(self.pref, current_window[1]), mode="w")

                last_time = format_time(current_window[0] - t0)
                print(
                    "Time t={}: {} reads and {} kbps.".format(
                        last_time, self.stats.nb_assigned_reads, round(self.stats.cumul_ln / 1000)
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
    """

    def __init__(self, rtbl):
        self.rtbl = rtbl
        self.phylogroups = sorted(self.rtbl.pgset.keys())
        self.summary = collections.OrderedDict()
        self.taxids = sorted(self.rtbl.pg.keys())

    def predict(self, stats):
        """Predict.
        """

        tbl = self.summary

        ## 1) GENERAL & CUMULATIVE STATS

        tbl['datetime'] = "NA"
        tbl['reads'] = stats.nb_assigned_reads + stats.nb_unassigned_reads
        tbl['bps'] = stats.cumul_ln
        tbl['matched bps'] = stats.cumul_h1

        ## 2) PG PREDICTION

        ## 2a) Find best and 2nd best PG
        ppgs = stats.pgs_by_weight()
        sorted_pgs = list(ppgs)
        pg1 = sorted_pgs[0]
        pg1_bm, pg1_w = ppgs[pg1]
        pg2 = sorted_pgs[1]
        pg2_bm, pg2_w = ppgs[pg2]

        ## 2b) Calculate PGS
        if pg1_w > 0:
            pgs = 2.0 * round(pg1_w / (pg1_w + pg2_w), 3) - 1
        else:
            pgs = 0.0

        ## 2c) Save values

        tbl['pgs'] = pgs
        tbl['pgs_ok'] = "pass" if pgs > 0.5 else "fail"
        tbl['pg1'] = pg1
        tbl['pg2'] = pg2
        tbl['pg1_bm'] = pg1_bm
        tbl['pg2_bm'] = pg2_bm
        tbl['pg1_w'] = round(pg1_w)
        tbl['pg2_w'] = round(pg2_w)

        ## 2c) Correct for missing data
        if pg1_w == 0:
            tbl['pg1'] = "NA"
            tbl['pg1_bm'] = "NA"
        if pg2_w == 0:
            tbl['pg2'] = "NA"
            tbl['pg2_bm'] = "NA"

        ## 3) ANTIBIOTIC RESISTANCE PREDICTION

        for ant in self.rtbl.ants:
            pres = stats.res_by_weight(pg1, ant)

            ##  3a) Calculate susceptibility score (sus) & correct for missing data
            try:
                # Identify S/R pivots
                s_bm = pres['S'][0]
                s_w = pres['S'][1]
                s_w_round = round(s_w)
                r_bm = pres['R'][0]
                r_w = pres['R'][1]
                r_w_round = round(r_w)
                # Calculate SUS
                if r_w + s_w > 0:
                    sus = round(s_w / (r_w + s_w), 3)
                else:
                    sus = 0.0
            except KeyError:
                # Some data were missing.
                if bm_cat == 'R':
                    # everything R
                    sus = 0.0
                    s_w, s_bm = "NA", "NA"
                elif bm_cat == 'S':
                    # everything S
                    sus = 1.0
                    r_w, r_bm = "NA", "NA"
                elif bm_cat == 'NA' and pg1_w == 0:
                    # not enough info yet (no match)
                    sus = 0.0
                    r_w, r_bm, s_w, s_bm = "NA", "NA", "NA", "NA"
                elif bm_cat == 'NA':
                    # ????
                    sus = 'NA'
                    r_w, r_bm, s_w, s_bm = "NA", "NA", "NA", "NA"

            ##  3b) Predict based on the collected info

            # best-match category
            if pg1_w > 0:
                bm_cat = self.rtbl.rcat[pg1_bm][ant]
            else:
                bm_cat = "NA"

            # prediction
            if sus > 0.6:
                pr_cat = "S"
            elif sus >= 0.5:
                pr_cat = "S!"
            else:
                pr_cat = "R"

            tbl[ant + "_sus"] = sus
            tbl[ant + "_pr_cat"] = pr_cat
            tbl[ant + "_bm_cat"] = bm_cat
            tbl[ant + "_r_bm"] = r_bm
            tbl[ant + "_s_bm"] = s_bm
            tbl[ant + "_r_w"] = r_w_round
            tbl[ant + "_s_w"] = s_w_round

            self.summary = tbl

    def print(self, timestamp):
        """Print.
        """

        global HEADER_PRINTED

        summary = self.summary

        if timestamp is None:
            current_dt, _, _ = str(datetime.datetime.now()).partition(".")
            summary['datetime'] = current_dt
        else:
            odt = datetime.datetime.fromtimestamp(timestamp)
            original_dt = odt.strftime('%Y-%m-%d %H:%M:%S')
            summary['datetime'] = original_dt

        keys = list(summary)
        values = [summary[k] for k in summary]
        if not HEADER_PRINTED:
            print(*keys, sep="\t")
            HEADER_PRINTED = True
        print(*format_floats(*values), sep="\t")
        sys.stdout.flush()


class Stats:
    """Statistics for RASE predictions.

    The attributes contain all necessary information that are necessary for
    predicting in RASE.

    Params:
        tree_fn (str): Filename of the Newick tree.

    Attributes:
        nb_assigned_reads (int): Number of processed reads.
        nb_nonprop_asgs (int): Number of processed alignments (before propagation).
        nb_asgs (int): Number of processed alignments (after propagation).

        cumul_h1 (float): Cumulative hit count.
        cumul_ln (float): Cumulative read length.

        stats_ct (dict): isolate -> number of processed reads
        stats_h1 (dict): isolate -> weighted h1
        stats_ln (dict): isolate -> weighted qlen
    """

    def __init__(self, tree_fn, rtbl):
        self._rtbl = rtbl
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
        self.cumul_h1 = 0
        self.cumul_ln = 0

        # statistics for individual isolates, "_unassigned_" for unassigned
        self.stats_ct = collections.defaultdict(lambda: 0.0)
        self.stats_h1 = collections.defaultdict(lambda: 0.0)
        self.stats_ln = collections.defaultdict(lambda: 0.0)

    def weight(self, isolate):
        """Get weight of an isolate.
        """
        # currently we define weight as stats_h1, but can be changed
        return self.stats_h1[isolate]

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

    def pgs_by_weight(self):
        """Sort phylogroups by weight.

        Returns:
            OrderedDict: pg -> (pivot_isolate, weight)
        """

        d = collections.defaultdict(lambda: [None, -1])

        for isolate in self._isolates:
            if isolate == FAKE_ISOLATE_UNASSIGNED:
                continue
            pg = self._rtbl.pg[isolate]
            val = self.weight(isolate)

            if val > d[pg][1]:
                d[pg] = isolate, val

        l = list(d.items())
        l.sort(key=lambda x: x[1][1], reverse=True)
        return collections.OrderedDict(l)

    def res_by_weight(self, pg, ant):
        """Return heaviest R/S isolates for a given antibiotic.

        Returns:
            dict: category -> (taxid, val)
        """

        d = collections.defaultdict(lambda: (None, -1))

        for isolate in self._rtbl.pgset[pg]:
            val = self.weight(isolate)
            cat = self._rtbl.rcat[isolate][ant].upper()

            if val > d[cat][1]:
                d[cat] = isolate, val

        return dict(d)

    def update_from_read(self, asgs):
        """Update statistics from assignments of a single read.

        Params:
            asgs (dict): Assignments.
        """

        assert len(asgs) > 0, "No assignments provided"

        asg0 = asgs[0]
        is_assigned = asg0["assigned"]

        # is read is assigned to at least 1 node?
        ln = asg0['ln']
        if is_assigned:
            h1 = asg0['h1']

            # 1) update cumulative stats
            self.nb_assigned_reads += 1
            self.nb_nonprop_asgs += len(asgs)

            self.cumul_h1 += h1
            self.cumul_ln += ln

            # 2) update stats for individual isolates
            nnames = [asg["rname"] for asg in asgs]
            bestmatching_isolates = self._descending_isolates(*nnames)
            l = len(bestmatching_isolates)
            self.nb_asgs += l
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
            self.cumul_ln += ln
            self.stats_ct[FAKE_ISOLATE_UNASSIGNED] += 1
            self.stats_ln[FAKE_ISOLATE_UNASSIGNED] += ln

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
                    *format_floats(self.stats_ct[isolate], digits=0),
                    *format_floats(
                        1.0 * self.stats_ct[isolate] / self.nb_assigned_reads if self.nb_assigned_reads != 0 else 0.0,
                        digits=3
                    ),
                    *format_floats(self.stats_ln[isolate], digits=0),
                    *format_floats(self.stats_ln[isolate] / self.cumul_ln if self.cumul_ln != 0 else 0.0, digits=3),
                    *format_floats(self.stats_h1[isolate], digits=0),
                    *format_floats(self.stats_h1[isolate] / self.cumul_h1 if self.cumul_h1 != 0 else 0.0, digits=3),
                ]
            )

        table.sort(key=lambda x: float(x[5]), reverse=True)

        for x in table:
            print(*format_floats(*x, digits=5), sep="\t", file=file)


class RaseMetadataTable:
    """RASE metadata table.

    This class loads data from a file with a description of individual isolates
    and preprocesses them so that they can be accessed via its internal
    structures.

    Attributes:
        pg: taxid -> phylogroup
        pgset: phylogroup -> set of taxids
        rcat: taxid -> antibiotic -> category
        weight: taxid -> weight
        ants: List of antibiotics.
    """

    def __init__(self, tsv):

        self.pg = {}
        self.pgset = collections.defaultdict(set)

        self.rcat = collections.defaultdict(dict)

        self.weight = {}

        with open(tsv, 'r') as f:
            tsv_reader = csv.DictReader(f, delimiter='\t')

            ants = filter(lambda x: x.find("_mic") != -1, tsv_reader.fieldnames)
            self.ants = list(map(lambda x: x.replace("_mic", ""), ants))

            for x in tsv_reader:
                taxid = x['taxid']

                pg = x['phylogroup']
                self.pg[taxid] = pg
                self.pgset[pg].add(taxid)

                for a in self.ants:
                    cat = x[a + "_cat"]
                    assert cat in set(["NA", "S", "R", "s", "r"])
                    self.rcat[taxid][a] = cat


class RaseBamReader:
    """Iterator over all assignments of individual reads in a BAM/RASE file.

    Assumes a non-empty BAM file.

    Params:
        bam_fn (str): BAM file name.
        out_bam_fn (str): Output BAM file name.

    Attributes:
        t1 (int): Timestamp of first read
    """

    def __init__(self, bam_fn, out_bam_fn):
        self.assignment_reader = SingleAssignmentReader(bam_fn, out_bam_fn)
        self._buffer = []
        self._finished = False
        try:
            self._load_assignment()
        except StopIteration:
            error("The provided SAM/BAM stream does not contain any alignments.")
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
        try:
            asg = next(self.assignment_reader)
        except StopIteration:
            raise StopIteration
        self._buffer.append(asg)

    def _extract_t1(self):
        self.t1 = timestamp_from_qname(self._buffer[0]["qname"])


class SingleAssignmentReader:
    """Iterator over individual assignments in BAM/RASE.

    Assumes that it is possible to infer read lengths (either
    from the ln tag, base sequence or cigars).

    Params:
        bam_fn (str): BAM file name.
        out_bam_fn (str): Output BAM file name.

    Attributes:
        samfile (pysam.AlignmentFile): PySAM file object.
    """

    def __init__(self, bam_fn, out_bam_fn):
        try:
            self.samfile = pysam.AlignmentFile(bam_fn, "rb")
        except ValueError as e:
            error("SAM/BAM stream from '{}' could not be read:".format(bam_fn), str(e))

        if out_bam_fn is not None:
            try:
                self.output_bamfile = pysam.AlignmentFile(out_bam_fn, "wb", header=self.samfile.header)
            except ValueError as e:
                error("Output BAM file '{}' could not be created:".format(output_bam_fn), str(e))
        else:
            self.output_bamfile = None

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

        #
        # Inferring length:
        # 1) if tag ln present, use it
        # 2) if seq present, use len (seq)
        # 3) infer length from cigar
        #

        # 1) acquire a new alignment from SAM
        try:
            alignment = next(self.alignments_iter)
            if self.output_bamfile is not None:
                self.output_bamfile.write(alignment)
        except StopIteration:
            raise StopIteration

        # 2) re-compute basic read-related variables if necessary
        self.qname = alignment.qname
        if self.qname != self.last_qname:
            try:
                self.read_ln = alignment.get_tag("ln")
            except KeyError:  # a ln tag is not present
                if alignment.seq is not None:
                    self.read_ln = len(alignment.seq)
                else:
                    self.read_ln = alignment.infer_read_length()
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
        'tree_fn',
        type=str,
        metavar='<tree.nw>',
        help='RASE database newick tree',
    )

    parser.add_argument(
        'metadata_fn',
        type=str,
        metavar='<db.tsv>',
        help='RASE database table with metadata',
    )

    parser.add_argument('bam_fn', type=str, metavar='<in.asgs.bam>', help='input RASE assignments (- for stdin)')

    parser.add_argument(
        'out_bam_fn',
        type=str,
        metavar='<out.asgs.bam>',
        help='output RASE assignments (to archive)',
        default=[],
        nargs='?',
    )

    parser.add_argument('-p', type=str, dest='pref', metavar='STR', help="output dir for samplings", default=None)

    parser.add_argument(
        '-t',
        dest='mode',
        choices=['clock', 'read'],
        help='time extraction mode',
        default='auto',
    )

    parser.add_argument(
        '-s',
        type=int,
        dest='delta',
        metavar='INT',
        help='sampling interval (seconds) [60]',
        default=60,
    )

    parser.add_argument(
        '-f',
        type=int,
        dest='first_read_delay',
        metavar='INT',
        help='delay of the first read (seconds) [60]',
        default=60,
    )

    args = parser.parse_args()

    if args.out_bam_fn:
        out_bam_fn = args.out_bam_fn
    else:
        out_bam_fn = None

    r = Runner(
        metadata_fn=args.metadata_fn,
        tree_fn=args.tree_fn,
        bam_fn=args.bam_fn,
        pref=args.pref,
        mode=args.mode,
        delta=args.delta,
        first_read_delay=args.first_read_delay,
        out_bam_fn=out_bam_fn,
    )

    try:
        r.run()
    except KeyboardInterrupt:
        error("Keyboard interrupt")


if __name__ == "__main__":
    main()
