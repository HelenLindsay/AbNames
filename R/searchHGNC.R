# searchHGNC ----

#'@title searchHGNC
#'
#'@description Search HGNC gene symbols, aliases and names for an exact match
#'to an value, assumed to be an antigen name or part of an antigen name.
#'
#'@param query_df A data.frame in long format for searching for exact matches to
#'HGNC names.  Must contain a column "value" containing potential gene symbols
#'for matching with HGNC symbols and descriptions and a column "ID" for grouping
#'results
#'@importFrom utils data
#'@importFrom dplyr left_join
#'@export
searchHGNC <- function(query_df){

    if (! all(c("ID", "value") %in% colnames(query_df))){
        stop("Query data frame must contain columns named ID and value")
    }

    official_nms <- c("HGNC_SYMBOL", "HGNC_NAME")

    utils::data("hgnc_long", envir = environment())
    utils::data("hgnc", envir = environment())

    res <- dplyr::left_join(query_df, hgnc_long) %>%
        dplyr::filter(! is.na(HGNC_ID)) %>%
        dplyr::group_by(ID) %>%

        # Count number of gene symbols per antibody ID
        dplyr::mutate(nsym_types = length(unique(symbol_type))) %>%

        # If there is more than one type of symbol, remove the previous symbols
        dplyr::filter(! (symbol_type == "prev_symbol" & nsym_types > 1)) %>%

        # If there is an exact match to the offical symbol, discard others
        dplyr::mutate(has_official = any(symbol_type %in% official_nms)) %>%
        dplyr::filter(has_official & symbol_type %in% official_nms |
                          ! has_official) %>%
        dplyr::select(-has_official) %>% #, -nsym_types) %>%

        # If there are matches to both symbol and name, keep symbol only
        dplyr::filter(! (any(symbol_type == "HGNC_SYMBOL") & !
                          symbol_type == "HGNC_SYMBOL")) %>%


        dplyr::mutate(nsym = length(unique(value)))

        # IF it is a multi-subunit protein, we expect all subunits to match






# Matches to greek word, subunit, upper_no_dash, missing matches to lower_no_dash?


    # If there are only matches to aliases/previous symbols, aggregate matches
    aliases <- res %>%
        dplyr::filter(! any(symbol_type %in% c("HGNC_SYMBOL", "HGNC_NAME"))) %>%

        dplyr::summarise(dplyr::summarise(HGNC_ALIAS = paste(unique(value), collapse = ", ")))


}

# CD77 should have matched previous symbol
# 4.1BB, 4.BBL should match dotToDash
# HLA.A.B.C?
# NKAT2 match to lower, Annexin V, Cadherin 11, Galectin-3
# HLA.A.B.C.


#ags %>%
#    dplyr::filter(! ID %in% res$ID) %>%
#    select(Antigen, Cat_Number) %>%
#    unique()
