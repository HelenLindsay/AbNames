# splitMerge ----

#' Apply a function to a subset of a data.frame
#'
#'@description Subset a data.frame according to a condition, apply a function
#' to the rows where the condition is TRUE, then rejoin with the rows where
#' condition is FALSE or NA.  A split-apply-combine where function is only
#' applied to a subset of rows.
#'
#' Filtering expression can be either quoted or unquoted, e.g.
#' complete.cases(x, y) (where x and y are column names) or
#' "complete.cases(x, y)".  To filter with multiple conditions use "&", e.g.
#' "x == 1 & y == 2"
#'
#' A use case for this function is when we want to perform a string processing
#' operation on only a subset of Antigens, for example, ignore isotype controls
#' or apply special formatting to only TCR antigens.
#'@param df A data.frame or tibble
#'@param ex character(1) An character expression for filtering df using
#'dplyr::filter, e.g. 'grepl("X", colname)'
#'@param f  A function to apply to the rows of df where ex is TRUE and that
#'returns a data.frame
#'@param verbose Should a warning be issued if extra rows are added after
#'applying f? (Default: TRUE)
#'@param ... Extra arguments for f
#'@importFrom rlang parse_expr enexpr is_string expr
#'@export
#'@returns df where function f has been applied only to the rows where
#'ex is TRUE
#'@examples
#' # Create a formatting function to selectively alter Antigen column of df
#' f <- function(df, fmt){
#'    df$Antigen <- sprintf(fmt, df$Antigen)
#'    return(df)
#' }
#'
#'df <- data.frame(Antigen = c("CD4","HLA-abc", "Siglec-5"))
#'
#'Apply formatting function f to entries of df$Antigen that don't match "CD"
#'splitMerge(df, "! grepl('CD', Antigen)", f, fmt = "%s TRIVIAL EXAMPLE")
splitMerge <- function(df, ex, f, verbose=TRUE, ...){

    # Switch depending on whether ex is a string or an expression
    enex <- rlang::enexpr(ex)

    if (rlang::is_string(enex)){
        ex <- rlang::parse_expr(ex) # Parse string into expression
    } else {
        ex <- enex
    }

    tmp <- .tempColName(df)
    original_nrow <- nrow(df)
    df <- dplyr::mutate(df, !!tmp := !! ex)

     # Remove negative cases
    df_not_ex <- dplyr::filter(df, ! (!!sym(tmp)) | is.na(!!sym(tmp)))

    # Filter for positive case:
    df <- dplyr::filter(df, !! rlang::sym(tmp))

    # Apply f to filtered df and re-join
    df <- f(df, ...)

    result <- dplyr::full_join(df, df_not_ex,
                               by=intersect(colnames(df),
                                            colnames(df_not_ex))) %>%
        dplyr::select(-all_of(tmp))

    if (isTRUE(verbose) & ! nrow(result) == original_nrow){
        warning("Rows were added when merging split data.frames")
    }

    return(result)
}
