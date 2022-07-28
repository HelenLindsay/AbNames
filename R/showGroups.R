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
#'@param interactive (default: TRUE) Should the function wait for a
#'command prompt to show the next group?
#'@importFrom dplyr group_rows
#'@export
showGroups <- function(df, i = 1, n = 1, max_rows = 50, interactive = TRUE){
    if (! "grouped_df" %in% class(df)){
        warning("Data.frame is not grouped")
    }

    stop_interactive <- FALSE
    msg <- "Enter\nn to print the next group, or\nq to quit"
    msg2 <- "Please enter either n (next) or q (quit)"
    gp_info <- "Group %s of %s: %s rows\n"

    row_idxs <- df %>% dplyr::group_rows()
    n_groups <- length(row_idxs)

    while (isFALSE(stop_interactive)){
        # Print out group number i
        cat(sprintf(gp_info, i, n_groups, length(row_idxs[[i]])))
        p_df <- .getGroups(df, i = i, n = n, row_idxs = row_idxs)
        print(as.data.frame(p_df[seq_len(min(max_rows, nrow(p_df))), ]))
        i <- i + 1

        if (! isTRUE(interactive)) break()

        if (i > n_groups){
            cat("\nNo more groups to show\n")
            break()
        }

        choice <- readline(msg)

        if (! choice %in% c("n", "q")) readline(msg2)
        if (choice == "q") stop_interactive <- TRUE

    }
}


# printGroupMatch ----
#
#' Select a group from a data.frame, format into a printable string
#'
#' Get the first group of a grouped data.frame matching a condition and format
#' output for printing an error message.
#'
#' @param df A data.frame or tibble
#' @param flt An (unquoted) expression for using with dplyr::filter
#' @importFrom rlang enquo
.printGroupMatch <-function(df, flt){
    multi_df <- df %>% dplyr::filter(!!rlang::enquo(flt))
    first_group <- .getGroups(multi_df)
    fg <- paste(capture.output(print(first_group)), collapse = "\n")
    return(fg)
}


# print_n ----
#
# print a data.frame n rows at a time
print_n <- function(df, n = 20){
    gp_info <- "Rows %s - %s (%s rows total)\n"
    msg <- "Enter\nn to print the next group, or\nq to quit"

    brks <- seq_len(nrow(df))[seq_len(nrow(df)) %% n == 0]

    # If last break equals number of rows, remove
    if (tail(brks, 1) == nrow(df)){
        brks <- head(brks, -1)
    }
    starts <- c(1, brks + 1)
    ends <- c(brks, nrow(df))

    # If last start equals number of rows, remove
    if (tail(starts, 1) == nrow(df)){
        starts <- head(starts, -1)
        ends <- head(ends, -1)
    }

    for (i in seq_along(starts)){
        print(i)
        cat(sprintf(gp_info, starts[i], ends[i], nrow(df)))
        print(df[starts[i]:ends[i],])
        choice <- readline(msg)
        if (! choice %in% c("n", "q")) readline(msg)
        if (choice == "q" ) break()
        if (i == length(starts)){
            message("no more groups to show")
        }
    }

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
#'@param row_idxs (Optional, default: NULL) Indices of rows
#'@importFrom dplyr group_rows
.getGroups <- function(df, i = 1, n = 1, row_idxs = NULL){
    if (is.null(row_idxs)){ row_idxs <- df %>% dplyr::group_rows() }

    if (i > length(row_idxs)){
        msg <- "Cannot print group %s, there are only %s groups in df"
        stop(sprintf(msg, i, length(row_idxs)))
    }

    row_idxs <- unlist(row_idxs[i:min((i + n - 1), length(row_idxs))])
    return(df[row_idxs, ])
}


# Find an antibody in a data.frame and return all aliases
# e.g. abAliases(df, "value == 'CD3'")
abAliases <- function(df, ex, by = "HGNC"){
    # Switch depending on whether ex is a string or an expression
    enex <- rlang::enexpr(ex)

    if (rlang::is_string(enex)){
        ex <- rlang::parse_expr(ex) # Parse string into expression
    } else {
        ex <- enex
    }

    res <- union_join(df, df %>%
                          dplyr::filter(!! ex) %>%
                          dplyr::select( {{ by }} ) )

    return(res)
}
