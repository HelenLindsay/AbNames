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


# .gsubNA ----
#
# Wrapper for gsub, returns NA (or option provided) if pattern was not matched
#
.gsubNA <- function(pattern, replacement, x, no_match = NA){
    res <- gsub(pattern, replacement, x)
    return(.noDups(res, x, no_match))
}


# .noDups ----
#
# Compares one vector to a reference, sets to NA if equal to the reference
#
.noDups <- function(x, reference, no_match = NA){
    x[x == reference] <- no_match
    return(x)
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


# .freducePartial ----
#' Create a list of partial functions and apply in sequence
#'
#' Apply a function f repeatedly to different columns of a data.frame.  Assumes
#' that ellipsis arguments contain (only) one vector entry and that function f
#' should be applied sequentially to each of these entries.  Does not work with
#' unquoted column names
#'
#' @param df A data.frame or tibble
#' @param f Function that returns a (mutated) data.frame
#' @param cls character(1) Name of entry in ... for iterating.
#' @param ... Extra arguments for f
#' @importFrom purrr partial
#' @importFrom magrittr freduce
.freducePartial <- function(df, f, cls, ...){
    dots <- list(...)

    stopifnot(cls %in% names(dots))

    col_vals <- dots[[cls]]
    dots[[cls]] <- NULL

    partials <- lapply(col_vals, function(x){
        args <- c(structure(c(f, x), names = c(".f", cls)), dots)
        do.call(purrr::partial, args)
    })

    result <- magrittr::freduce(df, partials)
    return(result)
}
