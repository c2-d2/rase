# RASE - Resistance-Associated Sequence Elements

<!-- vim-markdown-toc GFM -->

* [Installation](#installation)
* [File formats](#file-formats)
* [Related repositories](#related-repositories)
* [License](#license)
* [Contact](#contact)

<!-- vim-markdown-toc -->

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
  | `ks` | ks (proportion of matched *k*-mers) |
  | `ls` | lineage score |
  | `ls_pass` | lineage score interpretation, 1=_passing_ 0=_failing_ |
  | `ln`, `alt_ln` | predicted and alternative lineage |
  | `bm`, `bm_{prop}` | best-matching strain (nearest neighbor) and its properties |
  | `{ant}_ss` | susceptibility score for the antibiotic `{ant}` |
  | `{ant}_pred` | prediction (score interpretation): `S`=susceptible, `R`=non-susceptible, `S!` and `R!`=low confidence calls |
  | `{ant}_bm` | resistance information for the best match, format: `cat (mic)` |

* **Prediction output (snapshot).** Tab-separated text file with the following columns:

  | column | description |
  | --- | --- |
  | `taxid` | taxid of a database strain, `_unassigned_` for reads without any k-mer matches with the database |
  | `lineage` | lineage |
  | `weight` | weight (cumulative "number of *k*-mers matches divided by the the number of best matches") |
  | `weight_norm` | normalized `weight` |
  | `length` | cumulative "read length divided by number of matches" |
  | `count` | cumulative "count divided by number of matches" |


## Related repositories

* [RASE supplementary](http://github.com/c2-d2/rase-supplement). Supplementary Materials for the RASE paper, including figures and tables.
* [ProPhyle](http://prophyle.github.io). A highly accurate and resource-frugal DNA sequence classifier used by RASE.
* [Prophex](http://github.com/prophyle/prophex). A k-mer index based on the Burrows-Wheeler Transform, used by ProPhyle.


## License

[MIT](LICENSE).


## Contact

[Karel Brinda](https://scholar.harvard.edu/brinda) \<kbrinda@hsph.harvard.edu\>

