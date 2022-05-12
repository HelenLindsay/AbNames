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

    utils::data("hgnc_long", envir = environment())
    utils::data("hgnc", envir = environment())

    res <- dplyr::left_join(query_df, hgnc_long) %>%
        dplyr::filter(! is.na(HGNC_ID)) %>%
        dplyr::group_by(ID) %>%

        # Count number of gene symbols per antibody ID
        dplyr::mutate(nsym = length(unique(value)),
                      nsym_types = length(unique(symbol_type))) %>%

        # If there is more than one type of symbol, remove the previous symbols
        dplyr::filter(! (symbol_type == "prev_symbol" & nsym_types > 1))


    # If there are only matches to aliases/previous symbols, aggregate matches
    aliases <- res %>%
        dplyr::filter(! any(symbol_type %in% c("HGNC_SYMBOL", "HGNC_NAME"))) %>%

        dplyr::summarise(dplyr::summarise(HGNC_ALIAS = paste(unique(value), collapse = ", ")))




    # Select the results that match a HGNC symbol or name (e.g. not an alias)
    official <- res %>%
        dplyr::filter(symbol_type %in% c("HGNC_NAME", "HGNC_SYMBOL")) %>%
        dplyr::select(-nsym, -symbol_type, -value, -nsym_types) %>%
        unique()

    # Join official symbols with aliases





}
