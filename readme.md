# RASE - Resistance-Associated Sequence Elements

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Installation](#installation)
* [File formats](#file-formats)
* [Related repositories](#related-repositories)
* [License](#license)
* [Contact](#contact)

<!-- vim-markdown-toc -->

## Introduction

This repository contains the core RASE software. Pipeline for users is located in a [separate repository](https://github.com/c2-d2/rase-pipeline/).

## Installation

1. Setup an [environment](https://github.com/c2-d2/rase-pipeline/blob/master/environment.md)

2. Install RASE using PIP

   ```bash
   pip install --upgrade git+https://github.com/c2-d2/rase
   ```

## File formats

* **Prediction output (timeline).**  Tab-separated text file with the following columns:

  | column | description |
  | --- | --- |
  | `datetime` | datetime of the read |
  | `reads` | number of processed reads |
  | `kbps` | number of processed bps (thousands) |
  | `kkmers` | number of matched *k*-mers (thousands) |
  | `ks` | *k*-mer score (proportion of matched *k*-mers) |
  | `ls` | lineage score |
  | `ls_pass` | lineage score interpretation, 1=_passing_ 0=_failing_ |
  | `ln`, `alt_ln` | predicted and alternative lineage |
  | `bm`, `bm_{feature}` | best-matching strain (nearest neighbor) and its features |
  | `{ant}_ss` | susceptibility score for the antibiotic `{ant}` |
  | `{ant}_pred` | prediction (score interpretation): `S`=susceptible, `R`=non-susceptible, `S!` and `R!`=low confidence calls |
  | `{ant}_bm` | resistance information for the best match, format: `cat (mic)` |

  See an [example file](tests/predict.tsv).

* **Prediction output (snapshot).** Tab-separated text file with the following columns:

  | column | description |
  | --- | --- |
  | `taxid` | taxid of a database strain, `_unassigned_` for reads without any k-mer matches with the database |
  | `lineage` | lineage |
  | `weight` | weight (cumulative "number of *k*-mers matches divided by the the number of best matches") |
  | `weight_norm` | normalized `weight` |
  | `length` | cumulative "read length divided by number of matches" |
  | `count` | cumulative "count divided by number of matches" |

  See an [example file](tests/snapshot.tsv).


* **RASE DB metadata.** Tab-separated text file with the following columns:

  | column | description |
  | --- | --- |
  | `taxid` | taxid of a database strain |
  | `lineage` | lineage |
  | `order` | order for plotting |
  | `{feature}` | arbitrary features (e.g., `serotype` or `MLST`)|
  | `{ant}_mic` | original MIC string (e.g., `<0.03`) |
  | `{ant}_int` | extracted MIC interval (`0-0.03`) |
  | `{ant}_cat` | resistance category  (`S`/`R`/`s`/`r`) |

  See an [example file](tests/metadata.tsv). Metadata files are generated in dedicated repositories (see [RASE DB skeleton](https://github.com/c2-d2/rase-db-skeleton), [*S. pneumoniae* RASE DB](https://github.com/c2-d2/rase-db-spneumoniae-sparc), and [*N. gonorrhoeae* RASE DB](https://github.com/c2-d2/rase-db-ngonorrhoeae-gisp)).


## Related repositories

* [RASE supplementary](http://github.com/c2-d2/rase-supplement). Supplementary Materials for the RASE paper, including figures, tables, experiments, and links to other related repositories.
* [RASE DB skeleton](http://github.com/c2-d2/rase-db-skeleton). Skeleton for creating novel RASE databases.


## License

[MIT](LICENSE).


## Contact

[Karel Brinda](https://scholar.harvard.edu/brinda) \<kbrinda@hsph.harvard.edu\>

