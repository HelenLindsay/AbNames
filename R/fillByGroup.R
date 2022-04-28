# fillByGroup ----
#' Fill NAs in data.frame by grouping
#'
#'Wrapper of tidyr::fill with checks for NAs in grouping values and the option
#'to fill with the majority value if there is more than one value per group.
#'If rlang is available and an error is raised, an example group that causes an
#'error is returned
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
    original_nrows <- nrow(df)

    # Remove rows with NAs in group, merge after filling
    na_rows <- df %>%
        dplyr::filter(!complete.cases(across(all_of(groups))))

    # Group data frame, check if there are multiple values per group
    df <- df %>%
        dplyr::filter(complete.cases(across(all_of(groups)))) %>%
        dplyr::group_by(!!!syms(groups)) %>%
        dplyr::mutate(!!temp_names[1] :=
                          dplyr::n_distinct(!!!syms(fill), na.rm = TRUE))

    # Case: multiple possible values in a fill group and we should stop
    if (multiple == "stop" & any(df[, temp_names[1]] > 1)){
        msg <- "Some fill columns have multiple values per group"
        if (requireNamespace("rlang")){
            multi_df <- df %>% dplyr::filter(!!sym(temp_names[1]) > 1)
            first_group <- .firstGroups(multi_df)
            rlang::abort(message = msg, df = first_group)
        } else {
            stop(msg)
        }
    }

    # Case: only one value per group
    if (all(df[, temp_names[1]] == 1)){
        df <- tidyr::fill(df, !!!syms(fill), .direction = "updown") %>%
            dplyr::ungroup()
    }

    # Case: multiple values per group, select the most frequent one
    if (multiple == "majority"){

    }

    df <- full_join(df, na_rows)
    if (! nrow(df) == original_nrows){
        warning("Rows have been lost or gained when merging NA rows")
    }
    return(df)
}