# .stopIfColExists ----
#'
#' Warn if a column already exists in a data.frame
#'
#'@param df A data.frame or tibble
#'@param new_col character(n) Name of a column(s) to add to df
.stopIfColExists <- function(df, new_col){
    # Check if new_col(s) already exists in df, stop if so

    if (any(new_col %in% colnames(df))){
        cols_exist <- toString(intersect(new_col, colnames(df)))
        msg1 <- sprintf("Column(s) %s already exists in data.frame.\n",
                        cols_exist)
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

# .wordGrep ----
#
# Grep requiring search string is a word
#
# Wrapper for grepl requiring that search string is either surrounded by spaces
# or appears at the start or the end of the search string.  Does not cover full
# stops at the end of sentences
.wordGrep <- function(pattern, x){
    grep_cmd <- sprintf("(^|\\s)%s($|\\s|\\.\\s|\\.$)", pattern)
    grepl(grep_cmd, x)
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


# groupsWith ----
#
# Select rows from a data.frame containing values in from another
#
#
#@param df1 Filtered data.frame
#@param df2 Unfiltered data.frame
#@param col character(n) Name of column to use for selecting rows from df1
#'@importFrom dplyr pull
groupsWith <- function(df1, df2, col){
    # Only keep the columns in df2
    col_vals <- df1 %>%
        dplyr::pull(!!sym(col)) %>%
        stats::na.omit()

    df2 <- df2 %>%
        dplyr::filter(!!sym(col) %in% col_vals)
    return(df2)
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



# union_join ----
# Note: names in "by" won't work as this isn't actually a join
#'@title Select rows matching any column in another data.frame
#'
#'@description Select values from a data.frame df matching any
#'column from another data.frame or a selection of row indices.  If a second
#'data frame and a selection of rows is provided, values from df matching
#'any value in df2[rows, ] are returned.  If only rows indices are provided,
#'rows matching any value in df[rows, ] are returned.
#'
#'@param df A data.frame from which to select matching rows
#'@param df2 Optional, a second data
#'@param rows Row indices for subsetting, either df2 if present or df
#'@param by = columns to select from df2
#'@keywords internal
union_join <- function(df, df2 = NULL, rows = NULL, by = NULL){
    if (! is.null(df2) & ! is.null(rows)){
        message("Row selection will be made from df2")
    }
    tmp <- .tempColName(df)
    qdf <- df
    if (! is.null(df2)) { qdf <- df2 }
    if (! is.null(rows)) qdf <- qdf[rows, ]
    if (is.null(by)) by <- colnames(qdf)

    # If there are names, assume they use dplyr join syntax, i.e. names
    # refer to columns in df, values to columns in df2
    if (! is.null(names(by))){
        no_name <- names(by) == ""
        names(by)[no_name] <- by[no_name]

        # Check that names of "by" exist in df
        if (! all(names(by) %in% colnames(df))){
            stop("Not all columns in 'by' appear in df:",
                    toString(setdiff(names(by), colnames(df))))
        }

        # Check that values of by exist in qdf
        if (! all(by %in% colnames(qdf))){
            stop("Not all columns in 'by' appear in df2:",
                 toString(setdiff(by, colnames(qdf))))
        }

        # Rename columns in qdf to match "by", set by to equal names of by
        qdf <- qdf[, by]
        colnames(qdf) <- names(by)
        by <- names(by)
    }

    if (! all(by %in% colnames(df))){
        warning("Not all columns in 'by' appear in df:",
                toString(setdiff(by, colnames(df))))
    }

    # Extra brackets needed, see https://github.com/tidyverse/dplyr/issues/6194
    df %>% dplyr::filter((dplyr::if_any(.cols = by,
                                .fns = ~.x %in% qdf[[dplyr::cur_column()]])))
}


# rm_ambiguous ----
# note doesn't check temp col name
# gp - name of column to group by
# id - name of column to check for duplications
rm_ambiguous <- function(df, gp, id){
    df %>% dplyr::group_by({{ gp }}) %>%
        dplyr::mutate(n_genes = dplyr::n_distinct( {{ id }} )) %>%
        dplyr::filter(n_genes == 1) %>%
        dplyr::select(-n_genes) %>%
        dplyr::ungroup()
}
