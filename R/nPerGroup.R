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
#'@importFrom rlang syms
nPerGroup <- function(df, group, col, nm = "n_per_group"){
    # To do? possible to vectorise col with n temp names

    # CHECK IF NAMES EXIST IN DF
    tmp <- .tempColName(df, nm = nm)

    df <- splitMerge(df, complete.cases(!!!syms(group)),
                     .addNPerGroup, group = group, nm = tmp, cols = col)
    return(df)
}


# group must be unquoted, cols works quoted and unquoted
.addNPerGroup <- function(df, group, cols){

    df <- df %>%
        dplyr::group_by({{ group }}) %>%
        dplyr::mutate(dplyr::across({{ cols }},
                                    .fns = list(ndistinct =
                                        ~dplyr::n_distinct(.x, na.rm = TRUE)),
                                    .names = "n{.col}")
                      )
    return(df)
}



#
#x <- citeseq %>%
#    dplyr::group_by(Cat_Number) %>%
#    dplyr::mutate(nrrid = n_distinct(RRID, na.rm = TRUE),
#                  nclone = n_distinct(tolower(Clone), na.rm = TRUE),
#                  ntotalseq = n_distinct(TotalSeq_Cat, na.rm = TRUE),
#                  noligo = n_distinct(Oligo_ID, na.rm = TRUE))

