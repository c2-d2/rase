# RASE - Resistance-Associated Sequence Elements

<!-- vim-markdown-toc GFM -->

* [Installation](#installation)
  * [Computational environment](#computational-environment)
    * [Dependencies](#dependencies)
    * [Setting up an environment](#setting-up-an-environment)
* [File formats](#file-formats)
* [Related repositories](#related-repositories)
* [License](#license)
* [Contact](#contact)

<!-- vim-markdown-toc -->

## Installation

### Computational environment

The easiest way how to setup a computational environment for RASE is using
[Bioconda](https://bioconda.github.io/). This approach has been tested on
multiple Unix and OS X machines, including clusters and virtual machines.

#### Dependencies

* [Python 3](https://www.python.org/downloads/)
* [ProPhyle](http://prophyle.github.io)
* [ETE 3](http://etetoolkit.org/)
* [PySAM](https://github.com/pysam-developers/pysam)
* [GNU Make](https://www.gnu.org/software/make/) or equivalent
* [GNU parallel](https://www.gnu.org/software/parallel/)
* [Ghost Script](https://www.ghostscript.com/)
* [SnakeMake](https://snakemake.readthedocs.io)
* [SAMtools](http://www.htslib.org/)
* [R](https://www.r-project.org/)
* [R OptParse](https://cran.r-project.org/web/packages/optparse/)
* [GCC 4.8+](https://gcc.gnu.org/) or equivalent
* [zlib](https://zlib.net/)


#### Setting up an environment

* **BioConda environment.** We recommend to create a separate software
  environment (here called `raseenv`):

    ```bash
    conda create -n raseenv prophyle ete3 pysam snakemake-minimal samtools parallel r-optparse
    source activate raseenv
    ```

* **BioConda default environment.** Alternatively, the packages can also be
  installed directly into the default BioConda environment. Nevertheless, this
  is not always reliable since some of the RASE dependencies might collide with
  packages that were installed previously.

    ```bash
    conda install prophyle ete3 pysam snakemake samtools parallel r-optparse pandas
    ```

* **Manually.** All the dependencies can also be installed without BioConda. Many
  of these packages are distributed using standard package systems such as
  [APT](https://wiki.debian.org/Apt).

    ```bash
    apt-get install build-essential python3 python3-pandas zlib1g-dev r-base r-cran-optparse ghostscript
    ```

  All the Python packages (ProPhyle, PySAM, ETE 3, and Snakemake) can be
  installed using [PIP](https://pypi.org/project/pip/):

    ```bash
    pip3 install prophyle pysam ete3 snakemake pandas
    ```


## File formats

* **Prediction output (timeline).**  Tab-separated text file with the following columns:

  | column | description |
  | --- | --- |
  | `datetime` | datetime of the read |
  | `reads` | number of processed reads |
  | `kbps` | number of processed bps (thousands) |
  | `kkmers` | number of matched k-mers (thousands) |
  | `kmers_prop` | proportion of matched k-mers |
  | `ls` | lineage score |
  | `ls_pass` | lineage score interpretation, 1=_passing_ 0=_failing_ |
  | `ln`, `alt_ln` | predicted and alternative lineage, respectively |
  | `bm`, `bm_{prop}` | best-matching isolate (nearest neighbor) and its properties |
  | `{ant}_ssc` | susceptibility score for the antibiotic `{ant}` |
  | `{ant}_pred` | its interpretation: `S`=susceptible, `R`=non-susceptible, `S!` and `R!`=low confidence calls |
  | `{ant}_bm` | resistance information for the best match, format: `cat (mic)` |

* **Prediction output (snapshot).** Tab-separated text file with the following columns:

  | column | description |
  | --- | --- |
  | `taxid` | taxid of a database isolate, `_unassigned_` for reads without any k-mer matches with the database |
  | `lineage` | lineage |
  | `weight` | weight (cumulative "number of k-mer best matches divided by the number of matches") |
  | `weight_norm` | normalized `weight` |
  | `length` | cumulative "read length divided by number of matches" |
  | `length_norm` | normalized `ln` |


## Related repositories

* [RASE supplementary](http://github.com/c2-d2/rase-supplement). Supplementary Materials for the RASE paper, including figures and tables.
* [ProPhyle](http://prophyle.github.io). A highly accurate and resource-frugal DNA sequence classifier used by RASE.
* [Prophex](http://github.com/prophyle/prophex). A k-mer index based on the Burrows-Wheeler Transform, used by ProPhyle.


## License

[MIT](LICENSE).


## Contact

[Karel Brinda](https://scholar.harvard.edu/brinda) \<kbrinda@hsph.harvard.edu\>

