# searchHGNC ----

#'@title Search for matches to antibody names in HGNC protein coding genes
#'
#'@description Search HGNC gene symbols, aliases and names for an exact match
#'to an value, assumed to be an antigen name or part of an antigen name.
#'query_df must have columns named "ID", "name" and "value".
#'
#'@param query_df A data.frame in long format for searching for exact matches to
#'HGNC names.  Must contain a column "value" containing potential gene symbols
#'for matching with HGNC symbols and descriptions, a column "ID" for grouping
#'results, and a column "name" indicating the type of match (only used for
#'multi-subunit proteins and TCRs).
#'@param multisubunit character(n) Name(s) of entries in "name" column
#'corresponding to multi-subunit proteins.  Use NA if none exist.
#'@return A data.frame of matches from IDs in query_df to HGNC genes, including
#'a column "n_matches" giving the number of matches to different genes.
#'@importFrom utils data
#'@importFrom dplyr left_join
#'@export
searchHGNC <- function(query_df, multisubunit = c("TCR_long", "subunit")){

    if (! all(c("ID", "name", "value") %in% colnames(query_df))){
        stop("Query data frame must contain columns named ID and value")
    }

    official_nms <- c("HGNC_SYMBOL", "HGNC_NAME")

    utils::data("hgnc_long", envir = environment())
    utils::data("hgnc", envir = environment())

    res <- dplyr::left_join(query_df, hgnc_long) %>%
        dplyr::filter(! is.na(HGNC_ID)) %>%

        # If it is a multi-subunit protein, we expect all subunits to match
        .checkSubunitMatches(query_df) %>%
        dplyr::group_by(ID) %>%

        # Count number of gene symbols per antibody ID
        dplyr::mutate(nsym_types = length(unique(symbol_type))) %>%

        # If there is more than one type of symbol, remove the previous symbols
        dplyr::filter(! (symbol_type == "prev_symbol" & nsym_types > 1)) %>%

        # If there is an exact match to the offical symbol, discard others
        dplyr::mutate(has_official = any(symbol_type %in% official_nms)) %>%
        dplyr::filter(has_official & symbol_type %in% official_nms |
                          ! has_official) %>%
        dplyr::select(-has_official, -nsym_types) %>%

        # If there are matches to both symbol and name, keep symbol only
        dplyr::filter(! (any(symbol_type == "HGNC_SYMBOL") & !
                          symbol_type == "HGNC_SYMBOL")) %>%

        # If there are only matches to aliases/previous symbols, aggregate
        group_by(ID, HGNC_ID, symbol_type) %>%
        dplyr::mutate(across(c(name, value),
                         ~ifelse(symbol_type %in% c("HGNC_SYMBOL", "HGNC_NAME"),
                                 .x, paste(unique(.x), collapse = "|")))) %>%
        unique() %>%

        dplyr::group_by(ID) %>%
        dplyr::mutate(n_matches = length(unique(value))) %>%
        dplyr::ungroup()

    return(res)
}
