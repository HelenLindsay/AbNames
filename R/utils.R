# .warnIfColExists ----
#'
#' Warn if a column already exists in a data.frame
#'
#'@param df A data.frame or tibble
#'@param new_col character(1) Name of a column to add to df
.warnIfColExists <- function(df, new_col){
    # Check if new_col already exists in df, stop if so

    if (new_col %in% colnames(df)){
        msg1 <- sprintf("Column %s already exists in data.frame.\n", new_col)
        msg2 <- "Please supply a different value for new_col"
        stop(msg1, msg2)
    }
}

# .tempColName ----
#' Create a name for a temporary column in a data.frame
#'
#' Returns a name that does not already exist in columns of a data.frame
#'@param df A data.frame or tibble
#'@param n (default 1) How many temporary column names are needed?
#'@param nm (character(1), default "TEMP") Prefix for temporary column names
.tempColName <- function(df, n = 1, nm = "TEMP"){
    cn <- colnames(df)
    i <- 1
    temp_col <- nm
    res <- character()

    while(length(res) < n){
        if (! nm %in% cn) res <- c(res, nm)
        nm <- sprintf("%s%s", temp_col, i)
        i <- i + 1
    }

    return(res)
}


# firstGroups ----
#' Return the first n groups of a grouped data.frame
#'
#' Useful for debugging, e.g. when trying to fill by but groups contain
#' multiple values.  Group order is determined by dplyr::group_rows
#'
#'@param df A grouped data.frame
#'@param n (Default: 1) How many groups should be returned?
#'@importFrom dplyr group_rows
.firstGroups <- function(df, n = 1){
    row_idxs <- df %>% dplyr::group_rows()
    row_idxs <- unlist(row_idxs[1:min(n, length(row_idxs))])
    return(df[row_idxs, ])
}


# .gsubNA ----
#'
# Wrapper for gsub, returns NA (or option provide) if pattern was not matched
#'
.gsubNA <- function(pattern, replacement, x, no_match = NA){
    res <- gsub(pattern, replacement, x)
    res[res == x] <- no_match
    return(res)
}


# .printf ----
#'
#' Wrapper for sprintf, doesn't convert NA to character
#'
#'@param pattern A pattern to use with sprintf
#'@param ... Extra arguments for sprintf
.printf <- function(pattern, ...){
    dots <- list(...)
    res <- sprintf(pattern, ...)
    any_na <- apply(data.frame(dots), 1, function(y) any(is.na(y)))
    res[any_na] <- NA
    return(res)
}

# splitMerge ----

#' Apply a function to a subset of a data.frame
#'
#' Subset a data.frame according to a condition, apply a function to the rows
#' where the condition is TRUE, then rejoin with the rows where condition is
#' FALSE.  A split-apply-combine where function is only applied to a subset of
#' rows.
#'
#'@param df A data.frame or tibble
#'@param ex character(1) An character expression for filtering df using
#'dplyr::filter, e.g. 'grepl("X", colname)'
#'@param f  A function to apply to the rows where ex is TRUE and returns a
#'data.frame
#'@param ... Extra arguments for f
#'@return df where function f has been applied only to the rows where ex is TRUE
#'@importFrom rlang parse_expr
.splitMerge <- function(df, ex, f, ...){

    # This works if input is an unquoted expression
    #df_ex <- dplyr::filter(df, {{ ex}} )
    #print(df_ex)
    #
    #df_not_ex <- dplyr::filter(df, ! {{ ex }})
    #print(df_not_ex)

    # Remove negative cases
    not_ex <- rlang::parse_expr( paste("!", ex))
    df_not_ex <- dplyr::filter(df, !!(not_ex))
    print(nrow(df_not_ex))

    # Filter for positive case works:
    ex <- rlang::parse_expr(ex)
    print(ex)

    df <- dplyr::filter(df, !!ex)
    print(nrow(df))


    df <- f(df, ...)

    result <- dplyr::full_join(df, df_not_ex)
    if (! nrow(result) == nrow(df) + nrow(df_not_ex)){
        warning("Rows were added when merging split data.frames")
    }

    return(result)
}
