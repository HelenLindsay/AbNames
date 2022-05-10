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
#'@format A data frame with 90788 rows and 7 variables:
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


# totalseq_cocktails ----

#' Antibody and gene information for TotalSeq cocktails from BioLegend
#'
#' This table contains information downloaded from www.biolegend.com and
#' preprocessed into a consistent format.  The TotalSeq A, B, C and D
#' (Heme Oncology) cocktails are included.
#'
#'@format A data frame with 485 rows and 8 variables:
#' \describe{
#'     \item{Oligo_ID}{BioLegend identifier for the oligo Barcode_Sequence}
#'     \item{Antigen}{Name of the antigen (antibody name minus "anti-" prefix)}
#'     \item{Clone}{Name of the cell line that produced the antibody}
#'     \item{Ensembl_ID}{Ensembl gene identifier}
#'     \item{Gene_Symbol}{Ensembl gene symbol}
#'     \item{TotalSeq_Cat}{Identifier for the TotalSeq cocktail}
#'     \item{Barcode_Sequence}{Sequence of DNA oligo conjuated to the antibody}
#'     \item{Reactivity}{Species in which the antibody will bind, usually human}
#' }
#'
#'@source \url{"https://www.biolegend.com/"}
"totalseq_cocktails"
