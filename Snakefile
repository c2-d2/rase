import glob
import os
import sys


snakemake.shell.prefix("set -eo pipefail;")
snakemake.shell("(./scripts/test_environments.sh 2>&1) >/dev/null || ./scripts/test_environments.sh")


# 1) Detect indexes

def file_size(fn):
    return os.stat(fn).st_size

def smallest_file(fns):
    size=10**30
    sfn=None
    for fn in fns:
        rfn=os.path.realpath(fn)
        if file_size(rfn)<size:
            size=file_size(rfn)
            sfn=fn
    return sfn

def find_db_files(suffix, smallest_only=False):
    fns=glob.glob("database/*.{}".format(suffix))
    if smallest_only:
        fns=[smallest_file(fns)]
    return [os.path.basename(fn).replace("."+suffix,"") for fn in fns]

indexes_base_tsv=find_db_files("tsv")

indexes_base_tar=find_db_files("tar.gz")
indexes=sorted(list(set(indexes_base_tsv) & set(indexes_base_tar)))

indexes_base_tar_smallest=find_db_files("tar.gz", smallest_only=True)
smallest_index=list(set(indexes_base_tsv) & set(indexes_base_tar_smallest))

print("Indexes:", indexes)
if len(indexes)==0:
    print("!!!! ", file=sys.stderr)
    print("!!!! WARNING ", file=sys.stderr)
    print("!!!! ", file=sys.stderr)
    print("!!!! No indexes files provided", file=sys.stderr)
    print("!!!! ", file=sys.stderr)


# 2) Detect experiments

fastqs=glob.glob("reads/*.fq")
experiments=sorted([os.path.basename(x[:-3]) for x in fastqs])
smallest_experiment=[os.path.basename(smallest_file(fastqs)[:-3])]
print("Experiments:", experiments)
if len(experiments)==0:
    print("!!!! ", file=sys.stderr)
    print("!!!! WARNING ", file=sys.stderr)
    print("!!!! ", file=sys.stderr)
    print("!!!! No FASTQ files provided", file=sys.stderr)
    print("!!!! ", file=sys.stderr)


localrules: test_environments

rule all:
    input:
        [
            [
                [
                    "prediction/.{}.fq.complete".format(e),
                    "prediction/.{}__{}.bam.complete".format(e, i),
                    "prediction/.{}__{}.quantify.complete".format(e, i),
                    "prediction/.{}__{}.predict.complete".format(e, i),
                    "plots/{}__{}.timeline.pdf".format(e, i),
                ]
                for i in indexes
            ]
            for e in experiments
        ],


rule test:
    input:
        [
            [
                [
                    "prediction/.{}.fq.complete".format(e),
                    "prediction/.{}__{}.bam.complete".format(e, i),
                    "prediction/.{}__{}.quantify.complete".format(e, i),
                    "prediction/.{}__{}.predict.complete".format(e, i),
                    "plots/{}__{}.timeline.pdf".format(e, i),
                ]
                for i in smallest_index
            ]
            for e in smallest_experiment
        ],


rule preprocess_reads:
    input:
        fq="reads/{pref}.fq",
    output:
        t="prediction/.{pref}.fq.complete"
    params:
        fq="prediction/{pref}.fq"
    benchmark:
        "benchmarks/{pref}.readprep.log"
    shell:
        """
            ./scripts/minion_rename_reads.py {input.fq} | paste -d '\t' - - - - | sort | uniq | perl -pe 's/\t/\n/g' > "{params.fq}"
            touch "{output.t}"
        """


rule classify:
    priority: 60
    input:
        "database/.{index}.complete",
        "prediction/.{pref}.fq.complete",
    output:
        t="prediction/.{pref}__{index}.bam.complete",
    params:
        fq="prediction/{pref}.fq",
        bam="prediction/{pref}__{index}.bam",
        index="database/{index}"
    benchmark:
        "benchmarks/{pref}__{index}.classify.log"
    shell:
        """
            prophyle classify -P "{params.index}" -m h1 "{params.fq}" \
                | samtools view -b > "{params.bam}"

            touch "{output.t}"
        """


rule quantify:
    priority: 70
    input:
        "prediction/.{pref}__{index}.bam.complete",
        "database/.{index}.complete",
    output:
        t="prediction/.{pref}__{index}.quantify.complete"
    params:
        bam="prediction/{pref}__{index}.bam",
    benchmark:
        "benchmarks/{pref}__{index}.quantify.log"
    shell:
        """
            mkdir -p "prediction/{wildcards.pref}__{wildcards.index}"
            scripts/prophyle_quantify.py -p 'prediction/{wildcards.pref}__{wildcards.index}/' -i 60 "database/{wildcards.index}/tree.nw" "{params.bam}" /dev/null
            touch "{output.t}"
        """


rule predict:
    priority: 80
    input:
        "prediction/.{pref}__{index}.quantify.complete",
        "database/.{index}.complete",
    output:
        t="prediction/.{pref}__{index}.predict.complete",
    params:
        tsv1="prediction/{pref}__{index}.predict_without_flags.tsv",
        tsv2="prediction/{pref}__{index}.predict.tsv",
    benchmark:
        "benchmarks/{pref}__{index}.predict.log"
    shell:
        """
            scripts/rase_predict.py database/{wildcards.index}.tsv prediction/{wildcards.pref}__{wildcards.index}/*.tsv > "{params.tsv1}"
            scripts/prediction_add_flags.py  {params.tsv1} > {params.tsv2}
            rm {params.tsv1}
            touch "{output.t}"
        """


rule plot_timeline:
    priority: 90
    input:
        "prediction/.{pref}__{index}.predict.complete",
        "database/.{index}.complete",
    output:
        pdf="plots/{pref}__{index}.timeline.pdf",
    benchmark:
        "benchmarks/{pref}__{index}.plot.log"
    params:
        tsv="prediction/{pref}__{index}.predict.tsv",
        index="database/{index}",
        pref="{pref}"
    shell:
        """
            ./scripts/plot_selected_snapshots.py "database/{wildcards.index}.tsv" "prediction/{wildcards.pref}__{wildcards.index}" 1 5 -1 "plots/{wildcards.pref}__{wildcards.index}.snapshots."
            for x in plots/{wildcards.pref}__{wildcards.index}.snapshots.*; do
                ./scripts/reproducible_pdf.sh "$x";
            done
            ./scripts/plot_timeline.R "{params.tsv}" "{output.pdf}"
            ./scripts/reproducible_pdf.sh "{output.pdf}"
        """


"""
Todo: add possible checking of index consistency
"""
rule decompress:
    output:
        t="database/.{index}.complete"
    input:
        gz="database/{index}.tar.gz",
        tsv="database/{index}.tsv",
    benchmark:
        "benchmarks/decompress.{index}.log"
    shell:
        """
            mkdir -p database/{wildcards.index}
            prophyle decompress {input.gz} database
            touch "{output.t}"
        """
