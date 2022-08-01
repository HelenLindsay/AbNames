# fillByGroup ----
#' Fill NAs in data.frame by grouping
#'
#'Wrapper of tidyr::fill with checks for NAs in grouping values and the option
#'to fill with the majority value if there is more than one value per group.
#'If an error is raised, an example group that causes an error is returned
#'
#'@param df A data.frame or tibble with missing (NA) values to be filled
#'@param group (character(n)) Names of column(s) to group by
#'@param fill (character(n)) Name(s) of column(s) to fill
#'@param multiple (Default: "stop") How should multiple values in columns to
#'be filled be handled?  Either "stop" (raise an error), "mode" (select
#'the most common value) or "ignore" (entries with multiple possible modes are
#'set to NA).
#'@param method Either "only_na", of only missing entries should be filled,
#' or "all", if all entries should be replaced with their group mode
#'@importFrom tidyr fill
#'@importFrom stats complete.cases
#'@importFrom dplyr group_by ungroup n_distinct coalesce all_of
#'@importFrom utils capture.output
#'@export
fillByGroup <- function(df, group, fill, method = c("only_na", "all"),
                        multiple = c("stop", "mode", "ignore")){
    .stopIfColExists(df, sprintf("n%s", fill))

    multiple = match.arg(multiple)
    overwrite = match.arg(method) == "all"
    original_nrows <- nrow(df)

    df <- AbNames::splitMerge(df, complete.cases(!!!syms(group)),
                              .fillByGroup, group = group, fill = fill,
                              multiple = multiple, overwrite = overwrite)

    if (! nrow(df) == original_nrows){
        warning("Rows have been lost or gained when merging NA rows")
    }

    return(df)
}


# .fillByGroup ----
#
#'@keywords internal
#'@importFrom dplyr all_of if_any
.fillByGroup <- function(df, group, fill, multiple, overwrite){
    # Group data frame, check if there are multiple values per group
    cn <- colnames(df)
    df <- .addNPerGroup(df, group, fill)
    n_per_gp <- setdiff(colnames(df), cn)

    # Case: multiple possible values in a fill group and we should stop
    if (multiple == "stop" & any(df[, n_per_gp] > 1)){
        msg <- "Some fill columns have multiple values per group, e.g. \n"
        fg <- .printGroupMatch(df, if_any(all_of(n_per_gp), ~">"(.x, 1)))
        stop(msg, fg)
    }

    # Case: only one value per group
    if (all(df[, n_per_gp] <= 1)){
        df <- tidyr::fill(df, !!!syms(fill), .direction = "updown") %>%
            dplyr::ungroup()
    }

    # Case: multiple values per group, select the most frequent one
    if (multiple %in% c("mode", "ignore")){
        df <- .freducePartial(df, groupMode, cls = "cl", cl = fill, gp = group,
                              keep_first = multiple == "mode",
                              overwrite = overwrite)
    }

    df <- dplyr::select(df, -all_of(n_per_gp))
    return(df)
}


# groupMode ----
#
#' Find the most common value per group
#'
#' Given a grouped data.frame, count values per group and return a vector with
#' the most common value for each group.  If there are several equally common
#' values, the first will be chosen.
#'
#'@param df a grouped tibble
#'@param cl Name of column to find mode
#'@param new_cl Name of column to create.  If NA (default), col is modified
#'@param gp Name(s) of columns to group by
#'@param min_n (integer(1), default NA) Minimum number of occurrences of
#'majority value.  If provided, the majority value will be set to NA when it
#'occurs less than min_n times.
#'@param keep_first (logical(1), default TRUE) If there are multiple modes,
#'should the first (in order of the data.frame) be selected?  If FALSE,
#'entries are set to NA
#'@param overwrite (logical(1), default FALSE) If FALSE, new_cl will only
#'replace NAs.  If TRUE, new_cl will contain the mode.  Note that
#'keep_first = FALSE with overwrite = TRUE may lead to values being replaced
#'with NA.
#'@importFrom dplyr n
#'@importFrom rlang .data
groupMode <- function(df, cl, gp, new_cl = NA, min_n = NA,
                      keep_first = TRUE, overwrite = FALSE){

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

    if (! isTRUE(keep_first)){
        # Set tmp to NA if there are multiple valid modes
        x <- .tempColName(df, 1, "x")
        df <- df %>%
            dplyr::mutate(!!x :=
                          sum(.data[[n]] == max(.data[[n]])) == max(.data[[n]]),
                      !!tmp := ifelse(.data[[x]] == TRUE, .data[[tmp]], NA)) %>%
        dplyr::select(-.data[[x]])
    }

    if (isFALSE(overwrite)){
        # Fill with the majority value if current value is NA
        df <- df %>%
            dplyr::mutate(!!new_cl :=
                          dplyr::coalesce(.data[[cl]], .data[[tmp]]))
    } else {
        # Set new_cl to the actual mode
        df <- df %>%
            dplyr::mutate(!!new_cl := .data[[tmp]])
    }

    df <- df %>%
        dplyr::select(-.data[[tmp]], -.data[[n]]) %>%
        dplyr::ungroup()
    return(df)

}

