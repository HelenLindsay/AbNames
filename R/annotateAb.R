#' Annotate antigens with identifiers
#'
#' This is a wrapper for the default annotation pipeline using the functions
#' makeQueryTable and searchAliases.  It works by first transforming antigen
#' names using a series of regular expressions, then looking for exact matches
#' of these in the gene_aliases table, and finally removing entries where
#' multiple genes were matched.
#'
#' If using this function, we recommend you check the results as the default
#' matching functions may not suit every data set. In particular, check the
#' results for HLA, TCR or Ig antigens.
#'@param x A data.frame or tibble.  Must include a column named "Antigen"
#'containing the antigen names
#'@param id_cols (character(n), default NA) Name(s) of columns use to identify
#'each row in x.
#'@param control_col (character(1), default NA) Name of a logical column in x
#'indicating whether an antigen is an isotype control.  Isotype controls will
#'not be matched.
#'@author Helen Lindsay
#'@returns x, with additional identifier columns including the database in
#'which the matching entry was found.
#'@export
#'@examples
#'# This example shows searching for alternative ways of writing the
#'# same antigen.
#'df <- data.frame(Antigen = c("CD279", "CD279 (PD-1)", "CD279_PD1",
#'"CD279(PD-1)", "PD-1 (CD279)","PD1 (CD279)"))
#'annotateAb(df)
annotateAb <- function(x, id_cols = NA, control_col = NA){
    if (is.na(id_cols)){ id_cols <- c("Antigen", "Study") }
    id_cols <- intersect(id_cols, colnames(x))

    x <- AbNames::addID(x, id_cols = id_cols)

    if (! "Antigen" %in% colnames(x)){
        stop("x must contain a column named 'Antigen'")
    }

    query_t <- x %>% dplyr::select(all_of(
        na.omit(c("ID", "Antigen", control_col))))

    query_t <- AbNames::makeQueryTable(query_t, ab = "Antigen",
                                       control_col = control_col)
    if (! is.na(control_col)){
        query_t <- query_t %>%
            dplyr::select(-dplyr::all_of(control_col))
    }

    idf_cols <- c("ALT_ID", "HGNC_ID", "HGNC_SYMBOL", "ENSEMBL_ID",
                 "UNIPROT_ID", "ENTREZ_ID", "SOURCE")

    alias_results <- searchAliases(query_t) %>%
        dplyr::group_by(.data$ID) %>%

        # Select ID cols
        dplyr::select(dplyr::matches("ID|HGNC|SOURCE|^name$")) %>%
        unique() %>% # Collapse results with same ID from different queries

        # Collapse multi-subunit entries, convert "NA" to NA
        dplyr::summarise(dplyr::across(all_of(idf_cols),
                                       ~toString(unique(.x)))) %>%
        dplyr::mutate(dplyr::across(all_of(idf_cols), ~dplyr::na_if(.x, "NA")))

    x <- x %>%
        dplyr::left_join(alias_results, by = "ID") %>%
        unique()
    if (any(is.na(x[["ALT_ID"]]))){ x <- searchTotalseq(x) }

    x <- x %>%
        dplyr::ungroup() %>%
        dplyr::select(-"ID")

   return(x)
}

