# nPerGroup ----
#
#'Group a data.frame and count unique values
#'
#'Grouping columns containing NAs are not considered.
#'
#'@param df a data.frame or tibble
#'@param group character(n) Character vector of name(s) of grouping columns
#'@param col character(n) Name(s) of column(s) in which to count distinct values
#'@importFrom dplyr n_distinct
#'@importFrom stats complete.cases
#'@importFrom rlang syms
nPerGroup <- function(df, group, col){

    # Stop if names to be added already exist
    .stopIfColExists(df, sprintf("n%s", col))

    df <- splitMerge(df, complete.cases(!!!syms(group)),
                     .addNPerGroup, group = group, cols = col)
    return(df)
}



# group and cols must be character
.addNPerGroup <- function(df, group, cols){
    if (identical(group, cols)) {
        stop("Column to group must differ from column to count")
    }

    df <- df %>%
        dplyr::group_by(dplyr::across(all_of(group))) %>%
        dplyr::mutate(dplyr::across(all_of(cols),
                                    .fns = list(ndistinct =
                                        ~dplyr::n_distinct(.x, na.rm = TRUE)),
                                    .names = "n{.col}")
                      )
    return(df)
}


