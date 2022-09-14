
# AbNames
CITE-seq (Cellular Indexing of Transcriptome and Epitopes by sequencing, [Stoeckius et al 2017](https://www.nature.com/articles/nmeth.4380))
is a method for simultaneously quantifying RNA and surface protein expression. CITE-seq uses antibodies conjugated to sequencing oligos to 
bind to cell surface proteins. The conjugated oligos are then sequenced along with the RNA, and are referred to as Antibody Derived Tags (ADTs).
There is presently no consensus for how to report ADTs.  For example, they may be reported by the antibody name or the name of the protein the antibody bound.   

AbNames is an R package for matching antibody names to gene, protein, or other identifiers,
and for standarding names across data sets.  We are aiming to submit to Bioconductor for the next release.  For more information, please see the vignettes.

Feel free to contact me with questions suggestions, or contributions.

## Installation

You can install the development version of AbNames from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("HelenLindsay/AbNames")
```
