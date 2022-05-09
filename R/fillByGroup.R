# fillByGroup ----
#' Fill NAs in data.frame by grouping
#'
#'Wrapper of tidyr::fill with checks for NAs in grouping values and the option
#'to fill with the majority value if there is more than one value per group.
#'If rlang is available and an error is raised, an example group that causes an
#'error is returned
#'
#'@param df A data.frame or tibble with missing (NA) values to be filled
#'@param group (character(n)) Names of column(s) to group by
#'@param fill (character(n)) Name(s) of column(s) to fill
#'@param multiple (Default: "stop") How should multiple values in columns to
#'be filled be handled?  Either "stop" (raise an error) or "mode" (select
#'the most common value)
#'@importFrom tidyr fill
#'@importFrom stats complete.cases
#'@importFrom dplyr group_by ungroup n_distinct coalesce all_of
#'@importFrom utils capture.output
#'@export
fillByGroup <- function(df, group, fill, multiple = c("stop", "mode")){
    multiple = match.arg(multiple)
    tmp <- .tempColName(df, 1, "ndistinct")
    original_nrows <- nrow(df)

    df <- AbNames::splitMerge(df, complete.cases(!!!syms(group)),
                              .fillByGroup, group = group, tmp = tmp,
                              fill = fill, multiple = multiple)

    if (! nrow(df) == original_nrows){
        warning("Rows have been lost or gained when merging NA rows")
    }

    return(df)
}


.fillByGroup <- function(df, group, tmp, fill, multiple){
    # Group data frame, check if there are multiple values per group
    df <- .addNPerGroup(df, group, tmp, fill)

    # Case: multiple possible values in a fill group and we should stop
    if (multiple == "stop" & any(df[, tmp] > 1)){
        msg <- "Some fill columns have multiple values per group, e.g. \n"
        fg <- .printGroupMatch(df, !!sym(tmp) > 1)
        stop(msg, fg)
    }

    # Case: only one value per group
    if (all(df[, tmp] == 1)){
        df <- tidyr::fill(df, !!!syms(fill), .direction = "updown") %>%
            dplyr::ungroup()
    }

    # Case: multiple values per group, select the most frequent one
    if (multiple == "mode"){
        if (length(fill) > 1){
            stop("Filling majority value for multiple ",
                 "columns is not implemented")
        }



        df <- groupMode(df, cl = fill, gp = group)
    }

    df <- dplyr::select(df, -all_of(tmp))
    return(df)
}


# groupMode ----
#
#' Find the most common value per group
#'
#' Given a grouped data.frame, count values per group and add a column with the
#' most common value for each group.  If there are several equally common
#' values, the first will be chosen.
#'
#'@param df a grouped tibble
#'@param cl Name of column find mode
#'@param new_cl Name of column to create.  If NA (default), col is modified
#'@param gp Name(s) of columns to group by
#'@param min_n (integer(1), default NA) Minimum number of occurrences of
#'majority value.  If provided, the majority value will be set to NA when it
#'occurs less than min_n times.
#'@param n Name of count column.
#'@importFrom dplyr n
#'@importFrom rlang .data
groupMode <- function(df, cl, gp, new_cl = NA, min_n = NA, n = NA){
    n <- .tempColName(df, 1, "n")
    tmp <- .tempColName(df, 1)

    if (is.na(new_cl)) { new_cl <- cl }

    # Find the mode for column cl in group gp
    df <- df %>%
        dplyr::group_by(!!!syms(c(gp, cl))) %>%
        dplyr::mutate(!!n := n(),
                      !!n := ifelse(is.na(!!sym(cl)), -1, !!sym(n))) %>%
        dplyr::group_by(!!!syms(gp)) %>%
        dplyr::mutate(!!tmp := .data[[cl]][.data[[n]] == max(.data[[n]])][1])

    if (! is.na(min_n)){
        df <- df %>%
            dplyr::mutate(!!tmp := ifelse(max(!!sym(n)) < min_n,
                                          NA, !!sym(tmp)))
    }

    # Fill with the majority value if current value is NA
    df <- df %>%
        dplyr::mutate(!!new_cl :=
                          dplyr::coalesce(.data[[cl]], .data[[tmp]])) %>%
        dplyr::select(-.data[[tmp]], -.data[[n]]) %>%
        dplyr::ungroup()

    return(df)
}

