# showGroups ----
#'
#'Interactively print individual groups of a grouped data.frame
#'
#'Print first group of a grouped data.frame, then prompt user to either print
#'the next group or quit.   Group order is determined by dplyr::group_rows.
#'
#'@param df A *grouped* data.frame
#'@param i (integer(1), default 1) Index of the group to be printed.
#'@param n (integer(1), default 1) How many groups to show at once?
#'@param max_rows integer(1) Maximum number of rows to print (Default: 50).
#'@importFrom dplyr group_rows
#'@export
showGroups <- function(df, i = 1, n = 1, max_rows = 50){
    stop_interactive <- FALSE
    msg <- "Enter\nn to print the next group, or\nq to quit"
    msg2 <- "Please enter either n (next) or q (quit)"
    gp_info <- "Group %s: %s rows"

    while (isFALSE(stop_interactive)){
        # Print out group number i
        p_df <- .getGroups(df, i = i, n = n)
        print(as.data.frame(p_df[1:min(max_rows, nrow(p_df)), ]))
        i <- i + 1

        choice <- readline(msg)

        if (! choice %in% c("n", "q")) readline(msg2)
        if (choice == "q") stop_interactive <- TRUE
    }
}


# printGpMatch ----
#
#' Select a group from a data.frame, format into a printable string
#'
#' Get the first group of a grouped data.frame matching a condition and format
#' output for printing an error message.
.printGroupMatch <-function(df, flt){
    multi_df <- df %>% dplyr::filter(!!enquo(flt))
    first_group <- .getGroups(multi_df)
    fg <- paste(capture.output(print(first_group)), collapse = "\n")
    return(fg)
}


# getGroups ----
#' Return the first n groups of a grouped data.frame
#'
#' Useful for debugging, e.g. when trying to fill by but groups contain
#' multiple values.  Group order is determined by dplyr::group_rows
#'
#'@param df A grouped data.frame
#'@param i (integer(1), default: 1) The index of the first group to return
#'@param n (integer(1), default: 1) How many groups should be returned?
#'@importFrom dplyr group_rows
.getGroups <- function(df, i = 1, n = 1){
    row_idxs <- df %>% dplyr::group_rows()

    if (i > length(row_idxs)){
        msg <- "Cannot print group %s, there are only %s groups in df"
        stop(sprintf(msg, i, length(row_idxs)))
    }

    row_idxs <- unlist(row_idxs[i:min((i + n - 1), length(row_idxs))])
    return(df[row_idxs, ])
}
