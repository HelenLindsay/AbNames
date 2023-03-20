# getCommonName ----

# default matching columns are Antigen, Cat_Number, Clone, HGNC_SYMBOL
# cols = columns for grouping
# ab = column for standardising
# fill_col = column to add
# n_matched = name of column to add with number of matching results
# ... pass keep = TRUE for keeping grouping columns for debugging
# Be careful of exceptions, e.g. Cat_Number == "custom made"
getCommonName <- function(x, cols = NULL, ab = "Antigen",
                          fill_col = "Antigen_std", n_matched = "n_matched",
                          ignore = list(Cat_Number = "[Cc]ustom"),
                          verbose = TRUE, ...){

    # TO DO: test for "n_matched" in x? At the moment this function
    # is only called internally, n_matched cannot be NULL

    keep_cols <- c(colnames(x), fill_col)

    if (fill_col %in% colnames(x)){
        stop(sprintf("Column %s already exists in data.frame", fill_col))
    }

    if (is.null(cols)) {
        cols <- c("Antigen", "Cat_Number", "Clone", "ALT_ID")
    }

    # Remove section after brackets
    if (! fill_col %in% colnames(x)){
        x <- dplyr::mutate(x, !! fill_col := gsub(" ?\\(.*", "", !!sym(ab)))
    }

    # Group by any e.g. catalogue number or exact match to antigen
    tmp_grp <- .tempColName(x, nm = "group")

    if (isTRUE(verbose)){
        message(sprintf("Grouping data using columns:\n%s",
                toString(cols)))
    }

    x <- group_by_any(x, groups = cols, new_col = tmp_grp, ignore = ignore) %>%
        dplyr::mutate(!!n_matched := dplyr::n())

    # Fill with most common value
    x <- fillByGroup(x, group = tmp_grp, method = "all",
                     multiple = "mode", fill = fill_col)

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
