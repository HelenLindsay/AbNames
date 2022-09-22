# getAliases ----
# Convenience function to return aliases

# Given an antibody (or gene) name, return all aliases in the gene_aliases
# table
# by Find aliases by matching ALT_ID (default) or HGNC_ID
# return - table of aliases, or nothing if no aliases found
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
