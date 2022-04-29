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
.tempColName <- function(df, n = 1){
    cn <- colnames(df)
    temp_col <- "TEMP"
    i <- 1
    res <- character()

    while(length(res) < n){
        if (! temp_col %in% cn) res <- c(res, temp_col)
        temp_col <- sprintf("TEMP%s", i)
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


# gsubNA ----
#'
# Wrapper for gsub, returns NA if pattern was not matched
#'
.gsubNA <- function(pattern, replacement, x){
    res <- gsub(pattern, replacement, x)
    res[res == x] <- NA_character_
    return(res)
}


# .printf ----
#'
#' Wrapper for sprintf, doesn't convert NA to character
#'
.printf <- function(pattern, ...){
    dots <- list(...)
    res <- sprintf(pattern, ...)
    any_na <- apply(data.frame(dots), 1, function(y) any(is.na(y)))
    res[any_na] <- NA
    return(res)
}

