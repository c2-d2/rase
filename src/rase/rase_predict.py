#! /usr/bin/env python3
"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

# ./scripts/rase_predict.py ~/github/my/rase-predict/database/spneumoniae_sparc.k18/tree.nw ~/github/my/rase-predict/database/spneumoniae_sparc.k18.tsv  ~/github/my/rase-predict/prediction/sp10_norwich_P33.filtered__spneumoniae_sparc.k18.bam | tL

# todo:
# - implement option auto for autodetection of the time mode
"""
    Architecture:

        class Worker:
            Wrapping everything together.

        class Predict:
            Predicting from assignment statistics.

        class Stats:
            Statistics for RASE predictions.

        class RaseDbMetadata:
            RASE DB metadata table.

        class RaseBamReader:
            Iterator over all assignments of individual reads (1 read = 1 record) in a BAM/RASE file.

"""

import argparse
import collections
import csv
import datetime
import ete3
import itertools
import glob
import json
import os
import pandas
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


def timestamp_from_datetime(date, time):
    """Get a timestamp from a date and time.

    Args:
        date (str): Date ('YY-MM-DD').
        time (str): Time ('hh:mm:ss').

    Returns:
        timestamp (int): Unix timestamp.
    """
    date_time = date + " " + time
    b = datetime.datetime.strptime(date_time, "%Y-%m-%d %H:%M:%S")
    ts = int(b.timestamp())
    return ts


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


class Worker:
    """
        Wrapping everything together.
    """

    def __init__(self, metadata_fn, tree_fn, bam_fn, out_bam_fn, pref,
                 final_stats_fn, mode, delta, first_read_delay, ls_thres_pass,
                 ss_thres_shiconf, ss_thres_sr, ss_thres_rhiconf,
                 mbp_per_min, mimic_datetime):
        self.mode = mode
        self.metadata = RaseDbMetadata(metadata_fn)
        self.stats = Stats(tree_fn, self.metadata)
        self.predict = Predict(
            self.metadata,
            ls_thres_pass=ls_thres_pass,
            ss_thres_shiconf=ss_thres_shiconf,
            ss_thres_sr=ss_thres_sr,
            ss_thres_rhiconf=ss_thres_rhiconf)
        self.rase_bam_reader = RaseBamReader(bam_fn, out_bam_fn)
        self.pref = pref
        self.final_stats_fn = final_stats_fn
        self.delta = delta
        self.first_read_delay = first_read_delay
        self.mbp_per_min = mbp_per_min

        self.mimic_datetime = mimic_datetime

    def run(self):
        # 1) set the initial window:
        #     [t0, t0+delta), where
        #        t0=t1-first_read_delay
        #        t1=time_of_first_read
        if self.mode == "clock":
            t1 = current_timestamp()
        elif self.mode == "read":
            t1 = self.rase_bam_reader.t1
        elif self.mode == "mimic-ont":
            t1 = timestamp_from_datetime(*self.mimic_datetime.split())
        else:
            assert 1 == 2, "Unknown mode provided ({})".format(self.mode)
        t0 = t1 - self.first_read_delay
        current_window = [t0, t0 + self.delta]  # [x, y)

        if self.pref is not None:
            f = open(
                "{}/{}.tsv".format(self.pref, current_window[1]), mode="w")

        # 2) iterate through individual reads, and update and print statistics
        for read_stats in self.rase_bam_reader:
            assert len(read_stats) > 0
            if self.mode == "clock":
                read_timestamp = current_timestamp()
            elif self.mode == "read":
                read_timestamp = timestamp_from_qname(read_stats[0]["qname"])
            elif self.mode == "mimic-ont":
                read_timestamp = t1 + (self.stats.cumul_ln /
                                       (10**6)) / (self.mbp_per_min / 60.0)
            else:
                assert 1 == 2, "Unknown mode provided ({})".format(self.mode)

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
                    f = open(
                        "{}/{}.tsv".format(self.pref, current_window[1]),
                        mode="w")

                last_time = format_time(current_window[0] - t0)
                print(
                    "Time t={}: {:,} reads and {:,} kbps.".format(
                        last_time, self.stats.nb_assigned_reads,
                        round(self.stats.cumul_ln / 1000)),
                    file=sys.stderr)

            self.stats.update_from_read(read_stats)

        if self.pref is not None:
            self.stats.print(file=f)
            f.close()

        if self.final_stats_fn is not None:
            with open(self.final_stats_fn, mode="w") as f:
                stats.print(file=f)


class Predict:
    """Predicting from assignment statistics.

    Loads prediction statistics for a given time point (relative
    similarity to individual samples).

    Attributes:
        metadata: Metadata table.
        phylogroups: Sorted list of phylogroups.
        summary: Summary table for the output.
        ls_thres_pass: Threshold for lineage passing.
    """

    def __init__(self, metadata, ls_thres_pass, ss_thres_shiconf,
                 ss_thres_sr, ss_thres_rhiconf):
        self.metadata = metadata
        self.lineages = sorted(self.metadata.lineageset.keys())
        self.summary = collections.OrderedDict()
        self.ls_thres_pass = ls_thres_pass
        self.ss_thres_shiconf = ss_thres_shiconf
        self.ss_thres_sr = ss_thres_sr
        self.ss_thres_rhiconf = ss_thres_rhiconf

    def predict(self, stats):
        """Predict.
        """

        tbl = self.summary

        ## 1) GENERAL & CUMULATIVE STATS

        tbl['datetime'] = "NA"
        tbl['reads'] = stats.nb_assigned_reads + stats.nb_unassigned_reads
        tbl['kbps'] = round(stats.cumul_ln / 1000)
        tbl['kkmers'] = round(stats.cumul_h1 / 1000)
        if stats.cumul_ln != 0:
            tbl['ks'] = round(stats.cumul_h1 / stats.cumul_ln, 5)
        else:
            tbl['ks'] = 0.0

        ## 2) LINEAGE PREDICTION

        ## 2a) Find best and alternative lineage
        lineages_by_w = stats.lineages_by_weight()
        sorted_lineages = list(lineages_by_w)
        if stats.cumul_ln > 0:
            lineage1 = sorted_lineages[0]
            lineage1_bm, lineage1_w = lineages_by_w[lineage1]
            lineage2 = sorted_lineages[1]
            lineage2_bm, lineage2_w = lineages_by_w[lineage2]

            ## 2b) Calculate lineage score
            if lineage1_w > 0:
                ls = 2.0 * round(lineage1_w / (lineage1_w + lineage2_w), 3) - 1
            else:
                ls = 0.0

        else:
            lineage1, lineage1_bm, lineage1_w, lineage2, lineage2_bm, lineage2_w = 6 * [
                "NA",
            ]
            ls = 0

        ## 2c) Save values

        tbl['ls'] = ls
        tbl['ls_pass'] = "1" if ls >= self.ls_thres_pass else "0"
        tbl['ln'] = lineage1
        tbl['alt_ln'] = lineage2

        if lineage1_w != "NA" and lineage1_w > 0:
            tbl['bm'] = lineage1_bm
        else:
            tbl['bm'] = "NA"

        for x in self.metadata.additional_cols:
            if stats.cumul_ln > 0:
                tbl["bm_" + x] = self.metadata.additional_info[lineage1_bm][x]
            else:
                tbl["bm_" + x] = "NA"

        ## 3) ANTIBIOTIC RESISTANCE PREDICTION

        for ant in self.metadata.ants:

            ## 3a) Find best-match category
            if lineage1_w != "NA" and lineage1_w > 0:
                bm_cat = self.metadata.rcat[lineage1_bm][ant]
                bm_mic = self.metadata.rmic[lineage1_bm][ant]
            else:
                bm_cat = "NA"
                bm_mic = "NA"

            pres = stats.res_by_weight(lineage1, ant)

            ##  3b) Calculate susceptibility score (ss) & correct for missing data

            # Identify S/R pivots
            try:
                s_bm, s_w = pres['S']
                s_w_round = round(s_w)
            except KeyError:
                s_bm, s_w, s_w_round = "NA", "NA", "NA"

            try:
                r_bm, r_w = pres['R']
                r_w_round = round(r_w)
            except KeyError:
                r_bm, r_w, r_w_round = "NA", "NA", "NA"

            # Calculate susceptibility scores
            if s_w != "NA" and r_w != "NA":
                if r_w + s_w > 0:
                    ss = round(s_w / (r_w + s_w), 3)
                else:
                    ss = 0.0
            elif s_w == "NA" and r_w == "NA":
                # no data yet
                ss = 0.0
            elif s_w == "NA" and r_w != "NA":
                # everything R
                assert bm_cat.upper() == "R" or r_w == 0, (bm_cat, pres)
                ss = 0.0
            elif s_w != "NA" and r_w == "NA":
                # everything S
                assert bm_cat.upper() == "S" or s_w == 0, (bm_cat, pres)
                ss = 1.0

            ##  3c) Predict based on the collected info

            # prediction
            if ss > self.ss_thres_shiconf:
                pr_cat = "S"
            elif ss > self.ss_thres_sr:
                pr_cat = "S!"
            elif ss > self.ss_thres_rhiconf:
                pr_cat = "R!"
            else:
                pr_cat = "R"

            tbl[f"{ant}_ss"] = ss
            tbl[f"{ant}_pred"] = pr_cat
            tbl[f"{ant}_bm"] = f"{bm_cat} ({bm_mic})"
            #tbl[ant + "_r_bm"] = r_bm
            #tbl[ant + "_s_bm"] = s_bm
            #tbl[ant + "_r_w"] = r_w_round
            #tbl[ant + "_s_w"] = s_w_round

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

    The attributes contain all information necessary for
    predicting in RASE.

    Params:
        tree_fn (str): Filename of the Newick tree.
        metadata (str): RASE DB metadata.

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

    def __init__(self, tree_fn, metadata):
        self._metadata = metadata
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
            descending_leaves[root.name] = set(
                [isolate.name for isolate in root])
        return descending_leaves

    def _descending_isolates(self, *nnames):
        isolates = set()
        for nname in nnames:
            isolates |= self._descending_isolates_d[nname]
        return sorted(isolates)

    def lineages_by_weight(self):
        """Sort phylogroups by weight.

        Returns:
            OrderedDict: lineage -> (pivot_isolate, weight)
        """

        d = collections.defaultdict(lambda: [None, -1])

        for isolate in self._isolates:
            if isolate == FAKE_ISOLATE_UNASSIGNED:
                continue
            lineage = self._metadata.lineage[isolate]
            val = self.weight(isolate)

            if val > d[lineage][1]:
                d[lineage] = isolate, val

        l = list(d.items())
        l.sort(key=lambda x: x[1][1], reverse=True)
        return collections.OrderedDict(l)

    def res_by_weight(self, lineage, ant):
        """Return heaviest R/S isolates for a given antibiotic.

        Returns:
            dict: category -> (taxid, val)
        """

        d = collections.defaultdict(lambda: (None, -1))

        for isolate in self._metadata.lineageset[lineage]:
            val = self.weight(isolate)
            cat = self._metadata.rcat[isolate][ant].upper()

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
                warnings.warn(
                    "A single read shouldn't be reported as unassigned multiple times ({})"
                    .format(asgs))
            self.nb_unassigned_reads += 1
            self.cumul_ln += ln
            self.stats_ct[FAKE_ISOLATE_UNASSIGNED] += 1
            self.stats_ln[FAKE_ISOLATE_UNASSIGNED] += ln

    def print(self, file):
        """Print statistics to a file.

        Args:
            file (file): Output file.
        """
        print(
            "taxid",
            "lineage",
            "w",
            "w_norm",
            "l",
            "l_norm",
            sep="\t",
            file=file)
        table = []
        for isolate in self._isolates + [FAKE_ISOLATE_UNASSIGNED]:
            if isolate == FAKE_ISOLATE_UNASSIGNED:
                lineage = "NA"
            else:
                lineage = self._metadata.lineage[isolate]
            table.append([
                isolate,
                lineage,
                *format_floats(self.stats_h1[isolate], digits=0),
                *format_floats(
                    self.stats_h1[isolate] / self.cumul_h1
                    if self.cumul_h1 != 0 else 0.0,
                    digits=3),
                *format_floats(self.stats_ln[isolate], digits=0),
                *format_floats(
                    self.stats_ln[isolate] / self.cumul_ln
                    if self.cumul_ln != 0 else 0.0,
                    digits=3),
                #*format_floats(self.stats_ct[isolate], digits=0),
                #*format_floats(
                #    1.0 * self.stats_ct[isolate] / self.nb_assigned_reads
                #    if self.nb_assigned_reads != 0 else 0.0,
                #    digits=3),
            ])

        table.sort(key=lambda x: int(x[2]), reverse=True)

        for x in table:
            print(*format_floats(*x, digits=5), sep="\t", file=file)


class RaseDbMetadata:
    """RASE DB metadata table.

    Loads RASE DB metadata file (describing individual DB isolates)
    and preprocesses them so that they can be accessed via its internal
    structures.

    Attributes:
        lineage: taxid -> lineage
        lineageset: lineage -> set of taxids
        rcat: taxid -> antibiotic -> category
        rmic: taxid -> antibiotic -> original_mic
        weight: taxid -> weight
        additional_info: taxid -> key -> value
        ants: List of antibiotics.
    """

    def __init__(self, tsv):

        self.lineage = {}
        self.lineageset = collections.defaultdict(set)

        self.rcat = collections.defaultdict(dict)
        self.rmic = collections.defaultdict(dict)

        self.additional_info = collections.defaultdict(dict)

        self.weight = {}

        df = pandas.read_csv(
            tsv,
            delimiter='\t',
            na_values=[],
            keep_default_na=False,
            dtype=str)
        df = df.rename(columns={'phylogroup': 'lineage', 'pg': 'lineage'})  # backward compatibility

        # extract antibiotic abbrev from col names
        re_mic = re.compile(r'(\w+)_mic')
        self.ants = re_mic.findall(" ".join(df.columns))

        print(
            "Antibiotics in the RASE DB:",
            ", ".join(self.ants),
            file=sys.stderr)

        # check that all required columns are present
        basic_cols = ["taxid", "lineage", "order"] + list(
            map(lambda x: x + "_cat", self.ants)) + list(
                map(lambda x: x + "_int", self.ants)) + list(
                    map(lambda x: x + "_mic", self.ants))
        self.additional_cols = set(df.columns) - set(basic_cols)
        for x in basic_cols:
            assert x in df.columns, f"Column '{x}' is missing in '{tsv}'."

        # check that values are ok
        cats = set(itertools.chain(*[df[ant + "_cat"] for ant in self.ants]))
        cats_wrong = cats - set(["S", "R", "s", "r"])
        assert not cats_wrong, "Unknown resistance categories: {}".format(
            ", ".join(cats_wrong))

        df_dict = df.to_dict('index')
        for _, row in df_dict.items():
            taxid = row['taxid']
            lineage = row['lineage']
            self.lineage[taxid] = lineage
            self.lineageset[lineage].add(taxid)

            for x in self.additional_cols:
                self.additional_info[taxid][x] = row[x]

            for ant in self.ants:
                self.rcat[taxid][ant] = row[ant + "_cat"]
                self.rmic[taxid][ant] = row[ant + "_mic"]


class RaseBamReader:
    """Iterator over all assignments of individual reads (1 read = 1 record) in a BAM/RASE file.

    Assumes a non-empty BAM file.

    Params:
        bam_fn (str): BAM file name.
        out_bam_fn (str): Output BAM file name.

    Attributes:
        t1 (int): Timestamp of first read
    """

    def __init__(self, bam_fn, out_bam_fn):
        self.assignment_reader = _SingleAssignmentReader(bam_fn, out_bam_fn)
        self._buffer = []
        self._finished = False
        try:
            self._load_assignment()
        except StopIteration:
            error(
                "The provided SAM/BAM stream does not contain any alignments.")
        self._extract_t1()

    def __iter__(self):
        return self

    def __next__(self):
        """Retrieve set of assignments of the next read.
        """

        if self._finished:
            raise StopIteration

        while len(self._buffer) < 2 or self._buffer[-1][
                "qname"] == self._buffer[-2]["qname"]:
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


class _SingleAssignmentReader:
    """Iterator over individual assignments (1 assignment = 1 record) in BAM/RASE.

    Assumes that it is possible to infer read lengths (either
    from the ln tag, base sequence or cigars).

    Params:
        bam_fn (str): BAM file name.
        out_bam_fn (str): Output BAM file name.

    Attributes:
        samfile (pysam.AlignmentFile): PySAM file object.
    """

    def __init__(self, bam_fn, out_bam_fn):
        self.samfile = None
        self.output_bamfile = None

        try:
            self.samfile = pysam.AlignmentFile(bam_fn, "rb")
        except ValueError as e:
            error("SAM/BAM stream from '{}' could not be read:".format(bam_fn),
                  str(e))

        if out_bam_fn is not None:
            try:
                self.output_bamfile = pysam.AlignmentFile(
                    out_bam_fn, "wb", header=self.samfile.header)
            except ValueError as e:
                error(
                    "Output BAM file '{}' could not be created:".format(
                        output_bam_fn), str(e))

        self.alignments_iter = self.samfile.fetch(until_eof=True)

        self.qname = None
        self.last_qname = None
        self.read_ln = None

    # todo: verify that it gets invoked when the program finishes prematurly due to using head in the commandline
    def __del__(self):
        if self.samfile is not None:
            self.samfile.close()
        if self.output_bamfile is not None:
            self.output_bamfile.close()

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
        except StopIteration:
            raise StopIteration

        if self.output_bamfile is not None:
            self.output_bamfile.write(alignment)

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
        help='RASE tree',
    )

    parser.add_argument(
        'metadata_fn',
        type=str,
        metavar='<db.tsv>',
        help='RASE metadata table',
    )

    parser.add_argument(
        'bam_fn',
        type=str,
        metavar='<in.asgs.bam>',
        help='input RASE assignments (- for stdin)')

    parser.add_argument(
        'out_bam_fn',
        type=str,
        metavar='<out.asgs.bam>',
        help='output RASE matches',
        default=[],
        nargs='?',
    )

    parser.add_argument(
        '-p',
        dest='pref',
        metavar='STR',
        help="output dir for samplings",
        default=None)

    parser.add_argument(
        '-f',
        dest='final_stats_fn',
        metavar='STR',
        help="statistics for the last snapshot",
        default=None)

    parser.add_argument(
        '-t',
        dest='mode',
        choices=['clock', 'read', 'mimic-ont'],
        help='time extraction mode',
        default='clock',
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
        '-0',
        type=int,
        dest='first_read_delay',
        metavar='INT',
        help='delay of the first read (seconds) [60]',
        default=60,
    )

    parser.add_argument(
        '--ls-thres-pass',
        type=float,
        dest='ls_thres_pass',
        metavar='FLOAT',
        help='lineage score threshold [0.5]',
        default=0.5,
    )

    parser.add_argument(
        '--ss-thres-shiconf',
        type=float,
        dest='ss_thres_shiconf',
        metavar='FLOAT',
        help='high-confidence S threshold [0.6]',
        default=0.6,
    )

    parser.add_argument(
        '--ss-thres-sr',
        type=float,
        dest='ss_thres_sr',
        metavar='FLOAT',
        help='S/R threshold [0.5]',
        default=0.5,
    )

    parser.add_argument(
        '--ss-thres-rhiconf',
        type=float,
        dest='ss_thres_rhiconf',
        metavar='FLOAT',
        help='high-confidence R threshold [0.4]',
        default=0.4,
    )

    parser.add_argument(
        '--mbp-per-min',
        type=float,
        dest='mbp_per_min',
        metavar='FLOAT',
        help='mbps per minute (for mimicking ONT) [0.5]',
        default=0.5,
    )

    parser.add_argument(
        '--datetime',
        dest='mimic_datetime',
        metavar='STR',
        help='datetime (for mimicking ONT) [2018-01-01 00:00:00]',
        default="2018-01-01 00:00:00",
    )

    args = parser.parse_args()

    if args.out_bam_fn:
        out_bam_fn = args.out_bam_fn
    else:
        out_bam_fn = None

    r = Worker(
        metadata_fn=args.metadata_fn,
        tree_fn=args.tree_fn,
        bam_fn=args.bam_fn,
        pref=args.pref,
        final_stats_fn=args.final_stats_fn,
        mode=args.mode,
        delta=args.delta,
        first_read_delay=args.first_read_delay,
        out_bam_fn=out_bam_fn,
        ls_thres_pass=args.ls_thres_pass,
        ss_thres_shiconf=args.ss_thres_shiconf,
        ss_thres_sr=args.ss_thres_sr,
        ss_thres_rhiconf=args.ss_thres_rhiconf,
        mbp_per_min=args.mbp_per_min,
        mimic_datetime=args.mimic_datetime,
    )

    try:
        r.run()
    except KeyboardInterrupt:
        error("Keyboard interrupt")


if __name__ == "__main__":
    main()
