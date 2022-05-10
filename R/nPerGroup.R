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
#'@importFrom lang syms
nPerGroup <- function(df, group, col, tmp){

    df <- splitMerge(df, complete.cases(!!!syms(group)),
                     .addNPerGroup, group = group, tmp = tmp, col = col)
    return(df)
}


.addNPerGroup <- function(df, group, tmp, col){
    df <- df %>%
        dplyr::group_by(!!!syms(group)) %>%
        dplyr::mutate(!!tmp := dplyr::n_distinct(!!!syms(col), na.rm = TRUE))
    return(df)
}
