# getCommonName ----

# default matching columns are Antigen, Cat_Number, Clone, HGNC_SYMBOL)
# cols = columns for grouping
# ab = column for standardising
# new_col = column to add
# ... pass keep = TRUE for keeping grouping columns for debugging
# Be careful of exceptions, e.g. Cat_Number == "custom made"
getCommonName <- function(x, cols = NULL, ab = "Antigen",
                          new_col = "Antigen_std",
                          ignore = list(Cat_Number = "[Cc]ustom"), ...){

    keep_cols <- c(colnames(x), new_col)

    if (new_col %in% colnames(x)){
        stop(sprintf("Column %s already exists in data.frame", new_col))
    }

    if (is.null(cols)) {
        cols <- c("Antigen", "Cat_Number", "Clone", "HGNC_SYMBOL", "ENSEMBL_ID")
    }

    # Antibody formatting makes this function less generalisable
    #x <- dplyr::mutate(x, !! new_col := !!sym(ab))
    ## Remove sections in brackets, replace Greek symbols
    #x <- dplyr::mutate(x, !! new_col := gsub("^[Aa]nti-| \\(.*", "", !!sym(ab)),
    #                   !! new_col := replaceGreekSyms(!!sym(new_col)))

    # Group by any e.g. catalogue number or exact match to antigen
    tmp_grp <- .tempColName(x, nm = "group")
    x <- group_by_any(x, groups = cols, new_col = tmp_grp, ignore = ignore)

    # Fill with most common value
    x <- fillByGroup(x, group = tmp_grp, method = "all",
                     multiple = "mode", fill = new_col)

    return(x)
}


fillByAny <- function(x, cols, fill, ignore = NULL, multiple = "mode",
                      method = "all", ...){
    dots <- list(...)

    tmp_grp <- .tempColName(x, nm = "group")
    x <- group_by_any(x, groups = cols, new_col = tmp_grp, ignore = ignore)

    # Fill with most common value
    x <- fillByGroup(x, group = tmp_grp, method = method,
                     multiple = multiple, fill = fill)

    if (isTRUE(dots$keep)){
        return(x)
    }

    # Remove temporary column
    return(dplyr::select(x, -tmp_grp))

}
