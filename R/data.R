# HGNC ----
#' Gene names, symbols and IDs from HGNC
#'
#' A table of gene ids, symbols, aliases, previous aliases, and names
#' from the Human Genome Naming Consortium, corresponding to genome build
#' GRCh38.
#'
#' @format A data frame with 30123 rows and 10 variables:
#' \describe{
#'     \item{HGNC_ID}{HGNC gene IDs}
#'     \item{HGNC_SYMBOL}{HGNC gene symbol}
#'     \item{HGNC_NAME}{HGNC gene name (short description in words)}
#'     \item{alias_symbol}{comma separated list of current gene symbol aliases}
#'     \item{prev_symbol}{comma separated list of formerly used symbol aliases}
#'     \item{ENSEMBL_ID}{Ensembl gene ID corresponding to HGNC_ID}
#'     \item{Chromosome}{Chromosomal location of the gene}
#'     \item{alias_name}{Names for gene aliases, separated by "|"}
#'     \item{prev_name}{Names for former gene aliases, separated by "|"}
#'     \item{UNIPROT_IDS}{Uniprot ID(s) corresponding to HGNC_ID}
#' }
#' @source \url{"https://www.genenames.org/cgi-bin/genegroup/download-all"}
"hgnc"

# HGNC (long format) ----
#' Gene ids, names and aliases from HGNC, in long format
#'
#' Contains the same data as hgnc, but aliases, previous aliases and names have
#' been split to contain one entry per row, with an additional "symbol_type" column
#' giving the source of the symbol (e.g. HGNC_SYMBOL, alias_symbol).  Chromosome
#' data is not included
#'
#'#' @format A data frame with 90788 rows and 7 variables:
#' \describe{
#'     \item{HGNC_ID}{HGNC gene IDs}
#'     \item{ENSEMBL_ID}{Ensembl gene ID corresponding to HGNC_ID}
#'     \item{alias_name}{Names for gene aliases, separated by "|"}
#'     \item{prev_name}{Names for former gene aliases, separated by "|"}
#'     \item{UNIPROT_IDS}{Uniprot ID(s) corresponding to HGNC_ID}
#'     \item{symbol_type}{Source of the "value" column, e.g. "HGNC_SYMBOL"}
#'     \item{value}{A gene symbol or name}
#'
#' }
"hgnc_long"
