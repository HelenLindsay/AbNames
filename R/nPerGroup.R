# nPerGroup ----
#
#'Group a data.frame and count unique values
#'
#'Grouping columns containing NAs are not considered.
#'
#'@param df a data.frame or tibble
#'@param group character(n) Character vector of name(s) of grouping columns
#'@param col character(1) Name of column in which to count distinct values
#'@param nm character(1), default: "n_per_group". Name of column containing
#'counts of distinct values per group.  If nm is already a column in df, an
#'integer suffix will be added.
#'@importFrom dplyr n_distinct
#'@importFrom stats complete.cases
nPerGroup <- function(df, group, col, nm = "n_per_group"){
    # To do?
    # It would be possible to vectorise col with n temp names

    dfo <- df

    .warnIfColExists(df, nm)

    # Remove rows with NAs in group
    na_rows <- dplyr::filter(df, !complete.cases(!!!syms(group)))

    # Group data frame, check if there are multiple values per group
    df <- dplyr::filter(df, complete.cases(!!!syms(group)))

    add_n <- function(df, group, col, nm){
        tmp <- .tempColName(df, nm = nm)
        df <- df %>%
            dplyr::group_by(!!!syms(group)) %>%
            dplyr::mutate(!!tmp :=
                              dplyr::n_distinct(!!!syms(col), na.rm = TRUE))
        return(df)
    }

    x <- add_n(df, group, col, nm)
    y <- .splitMerge()
}


