# gene_aliases ----
#' @title Gene names, symbols and IDs from HGNC, Entrez and Ensembl
#'
#' @description
#' A table of gene ids, symbols, aliases, previous aliases, and names
#' from the Human Genome Naming Consortium, NCBI (Entrez) and Ensembl,
#' corresponding to genome build GRCh38.  Data is primarily based on the HGNC
#' gene groups and protein-coding genes tables.  Additional, unambiguous gene
#' name aliases from Ensembl (fetched via Bioconductor package biomaRt) and NCBI
#' (fetched via the NCBI ftp site and the Bioconductor package org.Hs.eg.db)
#' have been added.  Mappings between HGNC, Ensembl and NCBI (Entrez) IDs are
#' mostly based on HGNC, with some corrections of obsolete Ensembl IDs using
#' the Ensembl data.  Ambiguous aliases, i.e. aliases shared by more than on
#' gene, have been removed.  Some ambiguous aliases may be protein names, but
#' others are abbreviations with different meanings.
#'
#' Aliases, previous aliases and names have been split to contain one
#' entry per row, with an additional "symbol_type" column giving the source of
#' the symbol (e.g. HGNC_SYMBOL, alias_symbol). The source tables have been
#' filtered by BIOTYPE to remove pseudogenes, read-through genes, RNA genes,
#' mitochondrial genes and genes of unknown biotype.  Only genes located on
#' chromosomes were fetched from Ensembl, i.e. not haplotypes or patches, with
#' the result that some genes do not have Ensembl IDs.
#'
#' @format A data frame with 131533 rows and 10 variables:
#' \describe{
#'     \item{HGNC_ID}{HGNC gene IDs}
#'     \item{ENSEMBL_ID}{Ensembl gene ID, from HGNC}
#'     \item{UNIPROT_ID}{UNIPROT ID, from HGNC}
#'     \item{HGNC_SYMBOL}{HGNC gene symbol}
#'     \item{ENTREZ_ID}{ENTREZ (NCBI gene) ID, from HGNC}
#'     \item{BIOTYPE}{Type of gene, usually from HGNC}
#'     \item{symbol_type}{Source of the "value" column, e.g. "HGNC_SYMBOL",
#'     "HGNC_NAME"}
#'     \item{value}{A gene symbol, symbol alias, or name}
#'     \item{ALT_ID}{HGNC ID or other relevant ID for a specific protein
#'     modification, isoform or carbohydrate.  If no stable ID was found,
#'     Antigen / Clone combination is used.}
#'     \item{SOURCE}{Source of the data.}
#' }
#' @source \url{https://www.genenames.org/cgi-bin/genegroup/download-all}
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt}
"gene_aliases"


# HGNC ----
#' @title Gene names, symbols and IDs from HGNC, in long format
#'
#' @description
#' A table of gene ids, symbols, aliases, previous aliases, and names
#' from the Human Genome Naming Consortium, corresponding to genome build
#' GRCh38.  Data is from the HGNC gene groups and protein-coding genes tables.
#' Aliases, previous aliases and names have been split to contain one
#' entry per row, with an additional "symbol_type" column giving the source of
#' the symbol (e.g. HGNC_SYMBOL, alias_symbol).
#'
#' The HGNC table has been filtered by BIOTYPE to remove pseudogenes,
#' read-through genes, RNA genes, mitochondrial genes and genes of unknown
#' biotype.
#'
#' @format A data frame with 109313 rows and 9 variables:
#' \describe{
#'     \item{HGNC_ID}{HGNC gene IDs}
#'     \item{ENSEMBL_ID}{Ensembl gene ID, from HGNC}
#'     \item{ENTREZ_ID}{ENTREZ (NCBI gene) ID, from HGNC}
#'     \item{UNIPROT_ID}{UNIPROT ID, from HGNC}
#'     \item{BIOTYPE}{Type of gene}
#'     \item{HGNC_SYMBOL}{HGNC gene symbol}
#'     \item{symbol_type}{Source of the "value" column, e.g. "HGNC_SYMBOL",
#'     "HGNC_NAME"}
#'     \item{value}{A gene symbol, symbol alias, or name}
#'     \item{SOURCE}{Source of the data.  Here always "HGNC"}
#' }
#' @source \url{https://www.genenames.org/cgi-bin/genegroup/download-all}
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt}
"hgnc"


# totalseq ----
#' Antibody and gene information for TotalSeq antibodies from BioLegend
#'
#' This table contains information downloaded from www.biolegend.com and
#' preprocessed into a consistent format.  Discontinued antibodies
#' do not appear in this table, and non-human gene identifiers have been
#' discarded.
#'
#' Note that we do not consider the mapping between antibodies and gene
#' identifiers to be very reliable.  We have removed cases where the same
#' antibody clone was mapped to different gene identifiers, but have not
#' checked that the mapping agrees with the HGNC information.
#'
#'@format A data frame with 1466 rows and 9 variables:
#' \describe{
#'     \item{Cat_Number}{BioLegend catalogue number}
#'     \item{Oligo_ID}{BioLegend identifier for the oligo Barcode_Sequence}
#'     \item{Antigen}{Name of the antigen (antibody name minus "anti-" prefix)}
#'     \item{Clone}{Name of the cell line that produced the antibody}
#'     \item{TotalSeq_Cat}{TotalSeq sequencing chemistry.  A, B, C or D.}
#'     \item{Reactivity}{Species in which the antibody will bind, usually Human}
#'     \item{Cross_Reactivity}{Species in which the antibody shows,
#'     cross-reactivity, where known}
#'     \item{Barcode_Sequence}{Sequence of DNA oligo conjuated to the antibody}
#'     \item{ENSEMBL_ID}{Ensembl gene identifier}
#' }
#'
#'@source \url{"https://www.biolegend.com/"}
"totalseq"
