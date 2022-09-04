# makeQueryTable ----

#'Separate Antibody/Antigen names into component parts
#'
#'Given a data.frame containing Antigen/Antibody names, reformat
#'the names into possible gene names and return the data.frame in long format.
#'
#'@param df A data.frame or tibble.
#'@param ab (character(1), default: "Antigen") Name of column in df containing
#'antibody names
#'@param id (default: "ID") Name of column in df containing IDs for each row.
#'@param control_col (character(1), default: NA)  Optional name of a logical
#'column indicating whether an antibody is an isotype control.  If present,
#'controls will be removed to avoid spurious matches as they usually do not
#'react against human genes.
#'@param fun Optional custom function for formatting antigen names.  Must take
#' an argument "ab" giving the column name (as above) as its only argument.
#' Other arguments can be prefilled e.g. with purrr::partial.
#'@importFrom tidyr pivot_longer
#'@export
makeQueryTable <- function(df, ab = "Antigen", id = "ID",
                           control_col = NA, fun = NA){
    # Remove controls
    if (! is.na(control_col)) {
        df <- dplyr::filter(df, ! (!!sym(control_col)))
    }

    if (is.na(fun)) query_fun <- defaultQuery

    cn <- colnames(df)

    funs <- query_fun(ab = ab)
    df <- magrittr::freduce(df, funs)

    # Make a list of column names that are query terms
    query_cols <- c(ab, setdiff(colnames(df), cn)) # New columns + original ab
    query_cols <- setdiff(query_cols, id) # minus ID column

    # Convert to long format
    df <- qryToLong(df, query_cols)

    # Remove redundant entries
    df <- dplyr::group_by(df, !!sym(id), .data$value) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

    return(df)
}


# qryToLong ----
#
# Convert wide query table to long format
qryToLong <- function(df, query_cols){
    df <- df %>%
        tidyr::pivot_longer(cols = all_of(query_cols)) %>%
        dplyr::filter(!is.na(.data$value)) %>%
        unique()
    return(df)
}


# defaultQuery ----
#' Default list of functions for making an antibody query table
#'
#'@param ab (character(1), default "Antigen") Name of the column in df
#' containing Antigen/Antibody names
#'@export
defaultQuery <- function(ab = "Antigen"){

    # Default transformation sequence for making query table ----

    # Add a check that all columns to act on are either in the table or created
    # earlier in the pipeline?

    nc <- "Antigen_split"
    split_merge_str <- sprintf('grepl("TCR", %s)', nc)


    qry = list(purrr::partial(gsubAb, ab = !!ab), # Remove A/antis
               purrr::partial(gsubAb, ab = !!ab, pattern = "\\s[Rr]ecombinant"),
               purrr::partial(splitUnnest, ab = !!ab, new_col = !!nc),# Brackets
               purrr::partial(splitUnnest, ab = !!nc),  # Commas
               # / _ or . if at least 3 on each side and not TCR
               purrr::partial(splitUnnest, exclude = "TCR",
                        split = "(?<=[A-z0-9-]{3})[\\/_\\.](?=[A-z0-9-]{3,})"),
               purrr::partial(.reformatAb, ab = !!nc),
               purrr::partial(separateSubunits, ab = !!nc, new_col = "subunit"),
               purrr::partial(splitMerge, ex = !!split_merge_str, f = formatTCR,
                              tcr = "greek_letter", verbose = FALSE),
               purrr::partial(formatIg, ig = "greek_letter")
    )

    return(qry)
}


# addID ----
#'add ID column
#'
#'Adds an ID column to a data frame
#'
#'Pastes a group of columns together to form an ID column.  If pasted values do
#'not uniquely identify rows, adds a number to the end.
#'@param df A data.frame or tibble
#'@param id_cols (character(n)) Names of columns to paste to form ID column
#'@param new_col (character(1), default: "ID") Name of new ID column
#'@param warn (TRUE/FALSE, default: TRUE) If TRUE, warn if IDs are not unique
#'@param sep (Default: __) Delimiter to use for pasting columns to form ID
#'@return df with an extra ID column
#'
#'@importFrom dplyr mutate group_by all_of
#'@importFrom dplyr n syms row_number
#'@export
addID <- function(df, id_cols = c("Antigen", "Study"), new_col = "ID",
                  warn = TRUE, sep = "__"){

    # Check id_cols exist in df and new_col does not
    if (! all(id_cols %in% colnames(df))){ stop("All id_cols must be in df") }
    .stopIfColExists(df, new_col)

    df <- df %>%
        dplyr::mutate(!!new_col :=
                          do.call(paste, c(!!syms(id_cols), sep = sep))) %>%
        dplyr::group_by(!!sym(new_col)) %>%
        dplyr::mutate(n = n(),
                      i = row_number(),
                      !!new_col := ifelse(n == 1, !!sym(new_col),
                                     paste(!!sym(new_col), .data$i, sep = sep)))

    if (isTRUE(warn) & max(df[, "n"]) > 1){
        warning("ID columns do not uniquely identify rows, row numbers added.")
    }

    df <- df %>% dplyr::select(-.data$n, -.data$i)

    return(df)
}


# dplyr version ----
#addID2 <- function(df, id_cols = c(Antigen, Study), new_col = ID,
#                  warn = TRUE){
#
#    #.stopIfColExists(df, new_col)
#
#    df <- mutate(df, {{ new_col }} := do.call(paste,
#                                             c(across({{id_cols}}), sep = "-")))
#
#    if (isTRUE(warn)){
#        # Check that the ID column uniquely identifies rows
#        max_group_n <- group_by(df, {{ new_col }}) %>%
#            mutate(n = n()) %>%
#            pull(n) %>%
#            max()
#
#        if (max_group_n > 1){
#            warning("There are duplicate values in ID column")
#        }
#    }
#
#    return(df)
#}


# gsubAb ----
#'
#'Convenience function to remove a pattern from a column in a data.frame
#'
#'Remove a pattern from a column, either modifying in place or creating a new
#'column.  The default pattern removes the prefix "anti-" or "Anti-".
#'
#'@param df A data.frame or tibble
#'@param ab (character(1), default "Antigen) Name of the column to remove
#'prefixes from
#'@param pattern (character(1)) A regular expression for matching in column ab.
#'@param replacement (character(1)) Replacement value, default "" (i.e. remove)
#'@param exclude (Default: NA) - DOES NOTHING YET
#'@param restrict (Default: NA) - DOES NOTHING YET
#'@param new_col (character(1), default NA Name of the column to add to df.
#'If NA, column ab is modified
gsubAb <- function(df, ab = "Antigen", pattern = "[Aa]nti-", replacement = "",
                   exclude = NA, restrict = NA, new_col = NA){
    if (is.na(new_col)) new_col <- ab
    df <- dplyr::mutate(df, !!new_col := gsub(pattern, replacement, !!sym(ab)))

    # Restrict would have to be a filter expression, e.g. a particular study
    # Would need to use a temp column as in splitUnnest
    # Is this actually important?

    #@param exclude (character(1), default NA) a regex - do not split if ab
    #matches.
    #@param restrict (character(1), default NA) a regex replace only if restrict
    # matches in column ab

    return(df)
}


# upperSquish ----
#'
#' Convert values to uppercase and remove punctuation
#'
#' Convert a character vector to uppercase and remove dashes, spaces and dots.
#' Values are NA if the new value is identical to the original value and are
#' only converted if there the original value has 2 segments separated by,
#' punctuation, e.g. "IFN-g" would become "IFNG" but "IFN-g R alpha-chain" would
#' have new value NA.  Punctuation is not removed if it separates groups of
#' numbers, e.g. "CD3.1"
#'
#'@param ab (character(n)) A vector of strings to transform
#'@export
upperSquish <- function(ab){
    p1 <- "(^[A-z]+[0-9]?)[-\\. ]?([A-z]+)$" # e.g. CD3-A
    p2 <- "(^[A-z]+)[-\\. ]?([A-z0-9]+)$" # e.g. IL-2Rb
    x <- gsub(p1, "\\1\\2", ab)
    x <- toupper(gsub(p2, "\\1\\2", x))
    return(.noDups(x, ab))
}


# lowerNoDash ----
#'
#' Convert values to lowercase and replace - with space
#'
#' Convert a character vector to lowercase and replace - with spaces.
#' Values are NA if the new value is identical to the original value.
#'
#'@param ab (character(n)) A vector of strings to transform
#'@export
lowerNoDash <- function(ab){
    result <- gsub("-", " ", tolower(ab))
    return(.noDups(result, ab))
}


# dashNotDot ----
#'
#' Convert dashes (-) to dots(.)
#'
#' Convert a character vector to lowercase and replace "-" with "."
#' Values are NA if the new value is identical to the original value.
#'
#'@param ab (character(n)) A vector of strings to transform
#'@export
dashNotDot <- function(ab){
    return(.gsubNA("\\.", "-", toupper(ab)))
}


# .reformatAb ----
# Perform a default sequence of reformatting operations
.reformatAb <- function(df, ab = "Antigen"){
    result <- df %>%
        dplyr::mutate(greek_word = replaceGreekSyms(!!sym(ab), "sym2word"),
                      greek_letter = greekToLetter(!!sym(ab)),
                      upper_no_dash = upperSquish(.data$greek_letter),
                      lower_no_dash = lowerNoDash(.data$greek_letter),
                      dash_not_dot = dashNotDot(.data$greek_letter))
    return(result)
}



# TO DO: CDw for workshop CDs
