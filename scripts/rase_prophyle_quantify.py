#! /usr/bin/env python3

"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import collections
import os
import sys
import pysam
import json
import ete3


FAKE_ISOLATE_UNASSIGNED="_unassigned_"


def timestamp_from_qname(qname):
    return int(qname.partition("_")[0])


def format_time(seconds):
    minutes=seconds//60
    hours=minutes//60
    minutes-=60*hours
    return "{}h{}m".format(hours, minutes)


class Stats:
    """Statistics for an experiment.

    Params:
        tree_fn (str): Filename of the Newick tree.

    Attributes:
        nb_reads (int): Number of processed reads.
        nb_nonprop_asgs (int): Number of processed alignments (before propagation).
        nb_asgs (int): Number of processed alignments (after propagation).

        cumul_qlen (float): Cumulative weighted read length.
        cumul_h1 (float): Cumulative weighted hit count.
        cumul_c1 (float): Cumulative weighted coverage.

        stats_count (dict): isolate -> number of processed reads
        stats_qlen (dict): isolate -> weighted qlen
        stats_h1 (dict): isolate -> squared h1
        stats_c1 (dict): isolate -> squared c1
    """

    def __init__(self, tree_fn):
        self.tree=ete3.Tree(tree_fn, format=1)
        self.isolates=sorted([isolate.name for isolate in self.tree])
        self.descending_isolates=self.precompute_descendants(self.tree)

        # stats for assigned reads
        self.nb_assigned_reads=0
        self.nb_nonprop_asgs=0
        self.nb_asgs=0

        # stats for unassigned reads
        self.nb_unassigned_reads=0

        # cumulative statistics for assigned reads
        self.cumul_h1_pow1=0.0
        self.cumul_c1_pow1=0.0
        self.cumul_ln_pow1=0.0

        # statistics for individual isolates, "_unassigned_" for unassigned
        self.stats_h1_pow0=collections.defaultdict(lambda:0.0)
        self.stats_h1_pow1=collections.defaultdict(lambda:0.0)
        self.stats_c1_pow1=collections.defaultdict(lambda:0.0)
        self.stats_ln_pow1=collections.defaultdict(lambda:0.0)

    def precompute_descendants(self, tree):
        descending_leaves={}
        for root in list(tree.traverse())+[tree]:
            descending_leaves[root.name]=set([isolate.name for isolate in root])
        return descending_leaves

    def get_number_of_assigned_strains(self, asgs):
        l=0
        for asg in asgs:
            l+=len(self.descending_isolates[asg["rname"]])
        return l

    def update_from_one_read(self, asgs):
        """Update statistics from assignments of a single read.

        Params:
            asgs (dict): Assignments.
        """

        assert len(asgs)>0, "No assignments provided"

        is_assigned=asgs[0]["assigned"]

        if is_assigned:
            l=self.get_number_of_assigned_strains(asgs)

            self.nb_assigned_reads+=1
            self.nb_nonprop_asgs+=len(asgs)
            self.nb_asgs+=l

            for asg in asgs:
                nname=asg["rname"]
                descending_isolates=self.descending_isolates[nname]
                self.update_cumuls(h1=asg["h1"], c1=asg["c1"], ln=asg["ln"], weight=len(descending_isolates)/l)
                self.update_strain_stats(descending_isolates, h1=asg["h1"], c1=asg["c1"], ln=asg["ln"], l=l)

        else:
            assert len(asgs)==1, "A single read shouldn't be reported as unassigned mutliple times (error: {})".format(asgs)
            asg=asgs[0]
            self.nb_unassigned_reads+=1
            self.update_strain_stats([FAKE_ISOLATE_UNASSIGNED], h1=0, c1=0, ln=asg["ln"], l=1)

    def update_cumuls(self, h1, c1, ln, weight):
        self.cumul_h1_pow1+= h1 * weight
        self.cumul_c1_pow1+= c1 * weight
        self.cumul_ln_pow1+= ln * weight

    def update_strain_stats(self, isolates, h1, c1, ln, l):
        for isolate in isolates:
            self.stats_h1_pow0[isolate]+=1.0 / l
            self.stats_h1_pow1[isolate]+=1.0 * (h1 / l)
            self.stats_c1_pow1[isolate]+=1.0 * (c1 / l)
            self.stats_ln_pow1[isolate]+=1.0 * (ln / l)

    def print(self, file):
        """Print statistics.

        Args:
            file (file): Output file.
        """
        print("taxid", "count", "count_norm", "ln", "ln_norm", "h1", "h1_norm", "c1", "c1_norm", sep="\t", file=file)
        table=[]
        for isolate in self.isolates + [FAKE_ISOLATE_UNASSIGNED]:
            table.append(
                [
                    isolate,
                    self.stats_h1_pow0[isolate],
                    1.0*self.stats_h1_pow0[isolate]/self.nb_assigned_reads if self.nb_assigned_reads!=0 else 0,
                    self.stats_ln_pow1[isolate],
                    self.stats_ln_pow1[isolate]/self.cumul_ln_pow1 if self.cumul_ln_pow1!=0 else 0,
                    self.stats_h1_pow1[isolate],
                    self.stats_h1_pow1[isolate]/self.cumul_h1_pow1 if self.cumul_h1_pow1!=0 else 0,
                    self.stats_c1_pow1[isolate],
                    self.stats_c1_pow1[isolate]/self.cumul_c1_pow1 if self.cumul_c1_pow1!=0 else 0,
                ]
            )

        table.sort(key=lambda x: x[5], reverse=True)

        for x in table:
            print(*x, sep="\t", file=file)


class AssignmentBlockReader:
    """Iterator over blocks of ProPhyle assignments in a BAM file.

    Assumes a non-empty BAM file.

    Params:
        bam_fn (str): BAM file name.

    Attributes:
        t1 (int): Timestamp of first read
    """

    def __init__(self, bam_fn):
        self.assignment_reader=AssignmentReader(bam_fn)
        self._buffer=[]
        self._finished=False
        self._load_assignment()
        self._extract_t1()

    def __iter__(self):
        return self

    def __next__(self):
        """Get next block of assignments of the same read.
        """

        if self._finished:
            raise StopIteration

        while len(self._buffer)<2 or self._buffer[-1]["qname"]==self._buffer[-2]["qname"]:
            try:
                self._load_assignment()
            except StopIteration:
                self._finished=True
                buffer=self._buffer
                self._buffer=[]
                return buffer

        buffer=self._buffer[:-1]
        self._buffer=self._buffer[-1:]
        return buffer

    def _load_assignment(self):
        asg=next(self.assignment_reader)
        self._buffer.append(asg)

    def _extract_t1(self):
        self.t1=timestamp_from_qname(self._buffer[0]["qname"])


class AssignmentReader:
    """Iterator over individual ProPhyle assignments in a BAM file.

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

        self.qname=None
        self.last_qname=None
        self.read_ln=None

    def __iter__(self):
        return self

    def __next__(self):
        """Get next assignment.

        Returns:
            assignment (dict): A dict with the following keys: "rname", "qname", "qlen, "assigned", "h1", "c1".
        """

        # 1) acquire a new alignment from SAM
        try:
            alignment=next(self.alignments_iter)
        except StopIteration:
            raise StopIteration

        # 2) re-compute basic read-related variables if necessary
        self.qname=alignment.qname
        if self.qname!=self.last_qname:
            try:
                self.read_ln=alignment.get_tag("ln")
            except KeyError: # a ln tag is not present
                if alignment.seq!="*":
                    self.read_ln=len(alignment.seq)
                else:
                    #self.read_ln=read.infer_read_length()
                    self.read_ln=read.infer_read_length()
        self.last_qname=self.qname

        # 3) infer assignment-related variables
        if not alignment.is_unmapped:
            asg = {
                "rname": alignment.reference_name,
                "qname": alignment.qname,
                "assigned": True,
                "ln": self.read_ln,
                "h1": alignment.get_tag("h1"),
                "c1": alignment.get_tag("c1"),
            }
        else:
            asg = {
                "rname": None,
                "qname": alignment.qname,
                "assigned": False,
                "ln": self.read_ln,
                "h1": None,
                "c1": None,
            }

        return asg


def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument('tree',
            type=str,
            metavar='<tree.nw>',
        )

    parser.add_argument('bam',
            type=argparse.FileType('r'),
            metavar='<assignments.bam>',
        )

    parser.add_argument('tsv',
            type=str,
            metavar='stats.tsv',
        )

    parser.add_argument('-p',
            type=str,
            dest='pref',
            metavar='STR',
            required=True,
            help="Output dir for samplings"
        )

    parser.add_argument('-i',
            type=int,
            dest='delta',
            metavar='INT',
            help='sampling interval (in seconds) [300]',
            default=300,
        )

    parser.add_argument('-f',
            type=int,
            dest='first_read_delay',
            metavar='INT',
            help='delay of the first read [60]',
            default=60,
        )


    args = parser.parse_args()

    stats=Stats(args.tree)
    assignment_block_reader=AssignmentBlockReader(args.bam)

    # 1) set the initial window:
    #     [t0, t0+delta), where t0=time_of_first_read-first_read_delay
    t0=assignment_block_reader.t1-args.first_read_delay
    current_window=[t0, t0+args.delta] # [x, y)
    f=open("{}/{}.tsv".format(args.pref, current_window[1]), mode="w")

    # 2) iterate through individual reads, and update and print statistics
    for read_stats in assignment_block_reader:
        assert len(read_stats)>0
        read_timestamp=timestamp_from_qname(read_stats[0]["qname"])

        # do we have to shift the window?
        if read_timestamp >= current_window[1]:
            stats.print(file=f)
            f.close()
            while read_timestamp >= current_window[1]:
                current_window[0]+=args.delta
                current_window[1]+=args.delta
            f=open("{}/{}.tsv".format(args.pref, current_window[1]), mode="w")

            last_time=format_time(current_window[0]-t0)
            print("Time t={}: {} reads and {} non-propagated ({} propagated) assignments processed.".format(last_time, stats.nb_assigned_reads, stats.nb_nonprop_asgs, stats.nb_asgs), file=sys.stderr)

        stats.update_from_one_read(read_stats)

    stats.print(file=f)
    f.close()

    with open(args.tsv, mode="w") as f:
        stats.print(file=f)

if __name__ == "__main__":
    main()
