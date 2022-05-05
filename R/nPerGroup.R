# nPerGroup ----
#
#'Group a data.frame and count unique values
#'
#'Grouping columns containing NAs are not considered.
#'
#'@param df a data.frame or tibble
#'@param group character(n) Character vector of name(s) of grouping columns
#'@param col character(1) Name of column in which to count distinct values
#'@param na.rm
#'@importFrom dplyr n_distinct
#'@importFrom stats complete.cases
nPerGroup <- function(df, group, col){

    # Remove rows with NAs in group, merge after filling
    na_rows <- df %>%
        dplyr::filter(!complete.cases(!!!syms(group)))

    # Group data frame, check if there are multiple values per group
    df <- df %>%
        dplyr::filter(complete.cases(!!!syms(group))) %>%
        dplyr::group_by(!!!syms(group)) %>%
        dplyr::mutate(!!tmp := dplyr::n_distinct(!!!syms(fill), na.rm = TRUE))


}


