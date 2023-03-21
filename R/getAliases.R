#' getAliases ----
#'
#'Given an antibody name, return aliases
#'
#'Convenience function to return aliases all aliases in the gene_aliases table,
#'given a gene symbol or antibody name.  Aliases can be found by matching either
#'the ALT_ID (default) or HGNC_ID.  Matching by the ALT_ID means that isoforms
#'are not considered aliases.  For example CD45RA will not be an alias of
#'CD45RO.  Matching by HGNC_ID means that aliases corresponding to isoforms or
#'modifications of the same gene will all be returned.
#'
#'@param ab (character(1)) Name of the antibody/gene to match
#'@param by One of ALT_ID or HGNC_ID
#'@return A table of aliases, or nothing if no aliases are found
#'@author Helen Lindsay
#'@export
getAliases <- function(ab, by = c("ALT_ID", "HGNC_ID")){
    utils::data("gene_aliases", envir = environment())
    by <- rlang::sym(match.arg(by))
    ab_res <- gene_aliases %>%
        dplyr::filter(value == ab)
    if (nrow(ab_res) == 0){
        message(sprintf("%s not found in aliases table", ab))
        return()
    }
    ab_res <- ab_res %>% dplyr::pull(!!by)
    res <- gene_aliases %>%
        dplyr::filter(!!by %in% ab_res)
    return(res)
}


# abAliases ----
#' Find an antibody in a data.frame and return all aliases
#'
#' Filter a data frame by an expression (as expression or character) and
#' select all rows matching the value in the filtered column.
#' More general version of getAliases, as the filter function can use any
#' column e.g. abAliases(df, "value == 'CD3'").
#'@param df  A data.frame or tibble to filter
#'@param ex An filtering expression, as either a character or an expression
#'@param by Name of the column to use for selecting matching entries
#'(Default: "HGNC_ID")
#'@author Helen Lindsay
abAliases <- function(df, ex, by = "HGNC_ID"){
    # Switch depending on whether ex is a string or an expression
    enex <- rlang::enexpr(ex)

    if (rlang::is_string(enex)){
        ex <- rlang::parse_expr(ex) # Parse string into expression
    } else {
        ex <- enex
    }

    res <- filter_by_union(df, df %>%
                               dplyr::filter(!! ex) %>%
                               dplyr::select( {{ by }} ) )

    return(res)
}
