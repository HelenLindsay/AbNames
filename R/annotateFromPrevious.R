# fillByGroup ----
#' Fill NAs in data.frame by grouping
#'
#'Wrapper of tidyr::fill with checks for NAs in grouping values and the option
#'to fill with the majority value if there is more than one value per group.
#'
#'@param df A data.frame or tibble with missing (NA) values to be filled
#'@param groups (character(n)) Names of column(s) to group by
#'@param fill (character(n)) Name(s) of column(s) to fill
#'@param multiple (Default: "stop") How should multiple values in columns to
#'be filled be handled?  Either "stop" (raise an error) or "majority" (select
#'the most common value)
#'@importFrom tidyr fill
#'@importFrom stats complete.cases
#'@importFrom dplyr across all_of group_by n_distinct
#'@export
fillByGroup <- function(df, groups, fill, multiple = c("stop", "majority")){
    multiple = match.arg(multiple)
    temp_names <- .tempColName(df, 2)

    na_rows <- df %>%
        dplyr::filter(!complete.cases(across(all_of(groups))))
    df <- df %>%
        dplyr::filter(complete.cases(across(all_of(groups)))) %>%
        dplyr::group_by(!!!syms(groups)) %>%
        dplyr::mutate(!!temp_names[1] :=
                          dplyr::n_distinct(!!!syms(fill), na.rm = TRUE))

    if (multiple == "stop" & any(df[,temp_names[1]] > 1)){
        multi_df <- df %>% dplyr::filter(!!sym(temp_names[1]) > 1)
        gp_rows <- multi_df %>% group_rows()

        stop()
    }


}
