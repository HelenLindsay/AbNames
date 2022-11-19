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
