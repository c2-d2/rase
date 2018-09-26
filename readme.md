# RASE prediction pipeline


## Introduction

This repository contains the RASE prediction pipeline. The method uses lineage
calling to identify antibiotic resistant clones from nanopore reads. In our
[paper](https://www.biorxiv.org/content/early/2018/08/29/403204), we
demonstrate on the example of pneumococcus that, using this approach,
antibiotic resistance can be predicted within minutes from the start of
sequencing. Please, look at the paper for more information.

The RASE [Snakemake](https://snakemake.readthedocs.io/) workflow is specified
within a single [Snakefile](Snakefile). When executed, the pipeline first
detects the provided RASE database(s) (in the `database` directory) and
nanopore reads (in the `reads` directory), and generates all
`<db>`-`<experiment>` combinations. In practice, the most common scenario is
usally "1 db vs. many experiments". After the detection step, reads and
database are pre-processed: reads are sorted by time of sequencing and the
database gets uncompressed (i.e., the full internal ProPhyle k-mer index
restored).  Subsequently, nanopore reads from individual experiments are
compared to the database(s) using [ProPhyle](http://prophyle.github.io), and
isolate, phylogroup and resistance to individual antibiotics predicted - all
as a function of time.  Finally, the obtained time characteristics, as well as
rank plots for selected moments, are visualized using R.


## Quick example

The following example demonstrates the power of the streptococcal RASE with metagenomic reads.
The entire computation takes only 6m on a standard laptop (MacBook Pro). Note
that this is the experiment from Figure 3 (with human reads removed in-silico).

```bash
# add Bioconda channels, and create and activate an environment for RASE
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda env list | grep ^rase \
	|| conda create -y -n rase prophyle ete3 pysam snakemake samtools parallel r-optparse
source activate rase

# clone this repository
git clone https://github.com/c2-d2/rase-pipeline

# download the default database
(cd rase-pipeline/database \
	&& wget https://github.com/c2-d2/rase-db/releases/download/v01/spneumoniae_sparc.k18.tar.gz \
	&& wget https://github.com/c2-d2/rase-db/releases/download/v01/spneumoniae_sparc.k18.tsv)

# download minion reads from a metagenomic experiment
(cd rase-pipeline/reads \
	&& wget https://zenodo.org/record/1405173/files/sp10_norwich_P33.filtered.fq)

# run the pipeline
make -C rase-pipeline
```


## Installing of RASE

**Installing dependencies.** See [RASE computational
enviroment](https://github.com/c2-d2/rase/blob/master/environment.md).

**Cloning the RASE pipeline.**
You can clone this repository using git

```bash
git clone https://github.com/c2-d2/rase-pipeline
```

or download it as a [single .tar.gz
file](https://github.com/c2-d2/rase-pipeline/archive/master.tar.gz).

**Installing a RASE database.** A RASE database should be placed into the
directory `database`.  Every database consists of two files: a compressed
ProPhyle index (`.tar`) and a table with metadata for individual database
isolates (`.tsv`). To be properly detected, both of the files should have the
same base name.

The default RASE database (Streptococcus pneumoniae, k=18) can be downloaded
from [RASE DB releases](https://github.com/c2-d2/rase-db/releases). A custom
database can be constructed using scripts and workflows from the [RASE DB
repository](https://github.com/c2-d2/rase-db).

**Placing nanopore reads.** Nanopore reads should be placed into the `reads`
directory as a single `.fq` file per sequencing experiment. Please, check the
suffix: `.fastq` files are not currently detected. Also, the pipeline assumes
that the provided reads keep the original naming convention from ONT. Reads
that were used in the paper can be downloaded from
https://zenodo.org/record/1405173.


## Running RASE

**Running prediction locally.** The RASE pipeline can be executed by the Make
command (all steps, i.e., preprocessing, predicting, and plotting):

```bash
make
```

This will run the entire RASE pipeline, including plotting. For the testing and
debugging purposes, it might be sometimes useful to run RASE on a small
example. The following command will detect the smallest provided database and
experiment, and will run RASE only for this single combination:

```bash
make test
```

**Running prediction on a cluster.** When multiple experiments and multiple
databases are combined, it can be useful to parallelize the computation. In the
default setting, RASE supports the [Harvard O2
cluster](https://rc.hms.harvard.edu/#cluster) (Slurm-based), but the
[configuration file](cluster.json) can be easily adjusted. Submitting to a
cluster using
[Snakemake](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution)
can be done by `make cluster`.

**Exporting outputs.** Outputs of the pipeline can be exported to a single
`.tar` archive by `make export`.

**Cleaning.** In some situations, it might be useful to clean intermediate
files.

* `make clean` - cleans all prediction outputs and plots, but keeps ProPhyle
  outputs (the most time consuming step)
* `make cleanall` - removes previous and also ProPhyle outputs
* `git clean -fxd` - removes all files that are not part of the git repository
  (including the databases and reads)


**Others.** For the list of available subcommands, see the output of `make
help`.


## Files and directories

* `benchmarks/` - Snakemake benchmarks. As soon as any step of the pipeline
  gets finished, a log file with information about timing and memory
  consumption will appear here.
* `database/` - Source database files.
   - Each database consists of two files: `<db>.tar.gz` and `<db>.tsv`.
* `logs/` - Logs from job submission systems.
* `plots/` - Plotted figures.
   - `<experiment>__<db>.timeline.pdf` - Prediction as a function of time.
   - `<experiment>__<db>.snapshot.<time>.pdf` - Rank plot for selected times (1
	 minute, 5 minutes, last minute)
* `prediction/` - Prediction files.
   - `<experiment>__<db>.fq` - Nanopore reads after renaming and sorting by
	 timestamp.
   - `<experiment>__<db>.bam` - Reads mapped onto the phylogenetic tree using
	 ProPhyle.
   - `<experiment>__<db>/<timestamp>.tsv` - Cumulative weights calculated for
	 all isolates from the DB at the time; `h1` corresponds to the weights used
	 in the paper.
   - `<experiment>__<db>.predict.tsv` - Prediction timeline (each row
	 corresponds to one minute).
* `reads` - Nanopore reads (`<experiment>.fq`).
* `scripts` - RASE scripts.
* `tests` - Testing data for scripts.


## FAQs

> Can I run RASE on a laptop?

Yes, RASE is designed primarily for laptops; the memory requirements are low
(hundreds of MB) and the slowest step, read assignment, usually takes between
several minutes up to 2 hours, in dependence on the amount of sequencing data.

> Why do you then support submitting jobs to a cluster?

A cluster might be useful when many sequencing
experiments are to be processed at the same time, or a battery of databases needs to be
evaluated. In most of situations, we a laptop is sufficient.


## Related repositories

* [RASE DB](http://github.com/c2-d2/rase-db). Code for constructing RASE databases and the released databases.
* [RASE supplementary](http://github.com/c2-d2/rase-supplement). Supplementary Materials for the RASE paper, including figures and tables.
* [ProPhyle](http://prophyle.github.io). A highly accurate and resource-frugal DNA sequence classifier used by RASE.
* [Prophex](http://github.com/prophyle/prophex). A k-mer index based on the Burrows-Wheeler Transform, used by ProPhyle.

## Citing RASE

Karel Brinda, Alanna Callendrello, Lauren Cowley, Themoula Charalampous, Robyn
S Lee, Derek R MacFadden, Gregory Kucherov, Justin O'Grady, Michael Baym,
William P Hanage. **Lineage calling can identify antibiotic resistant clones
within minutes.**
bioRxiv, 2018.
doi:[10.1101/403204](https://doi.org/10.1101/403204)


## License

[MIT](LICENSE).


## Contact

[Karel Brinda](https://scholar.harvard.edu/brinda) \<kbrinda@hsph.harvard.edu\>

