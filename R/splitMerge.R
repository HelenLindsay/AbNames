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
splitMerge <- function(df, ex, f, ...){

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
    df <- dplyr::filter(df, !!ex)

    # Apply f to filtered df and re-join
    df <- f(df, ...)
    result <- dplyr::full_join(df, df_not_ex)

    if (! nrow(result) == nrow(df) + nrow(df_not_ex)){
        warning("Rows were added when merging split data.frames")
    }

    return(result)
}
