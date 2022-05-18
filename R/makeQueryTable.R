# makeQueryTable ----
#@param funs A vector of string formatting functions to apply to df.  Each
#should take a data.frame as a first argument and return a data.frame or
#tibble, and optional arguments should be filled in e.g. with purrr::partial.
#@param query_cols (character(n), default: NULL) Name(s) of columns containing
#potential names to search.  If not provided, any columns that are added to df
#after applying funs are used as query columns.
#
#'Separate Antibody/Antigen names into component parts
#'
#'Given a data.frame containing Antigen/Antibody names, reformat
#'the names into possible gene names and return the data.frame in long format.
#'
#'@param df A data.frame or tibble
#'@param ab (character(1), default: "Antigen") Name of column in df containing
#'antibody names
#'@param control_col (character(1), default: NA)  Optional name of a logical
#'column indicating whether an antibody is an isotype control.  If present,
#'controls will be removed to avoid spurious matches as they usually do not
#'react against human genes.
#'@importFrom tidyr pivot_longer
#'@export
makeQueryTable <- function(df, ab = "Antigen", control_col = NA){
    # Remove controls
    if (! is.na(control_col)) { df <- dplyr::filter(df, ! (!!sym(control_cl))) }

    cn <- colnames(df)

    funs <- defaultQuery(ab = ab)
    df <- magrittr::freduce(df, funs)

    # Make a list of column names that are query terms
    query_cols <- c(ab, setdiff(colnames(df), cn)) # New columns + original ab
    query_cols <- setdiff(query_cols, "ID") # minus ID column

    # Convert to long format
    df <- qryToLong(df, query_cols)

    return(df)
}


# qryToLong ----
#
# Convert wide query table to long format
qryToLong <- function(df, query_cols){
    df <- df %>%
        tidyr::pivot_longer(cols = all_of(query_cols)) %>%
        dplyr::filter(!is.na(value)) %>%
        unique()
    return(df)
}


# defaultQuery ----
#' Default list of functions for making an antibody query table
#'
#'@param ab (character(1), default "Antigen") Name of the column in df
#' containing Antigen/Antibody names
#'@param id_cols (character(n)) Names of columns to paste to form ID column
#'@export
defaultQuery <- function(ab = "Antigen"){

    # Default transformation sequence for making query table ----

    # Check that all columns to act on are either in the table or created
    # earlier in the pipeline?
    split_merge_str <- sprintf('grepl("TCR", %s)', ab)

    qry = list(purrr::partial(gsubAb, ab = !!ab), # Remove A/antis
               purrr::partial(gsubAb, ab = !!ab, pattern = "\\s[Rr]ecombinant"),
               purrr::partial(splitUnnest, ab = !!ab), # Brackets
               purrr::partial(splitUnnest, ab = !!ab, split = ", "),  # Commas
               # / _ or . if at least 3 on each side and not TCR
               purrr::partial(splitUnnest, exclude = "TCR",
                        split = "(?<=[A-z0-9-]{3})[\\/_\\.](?=[A-z0-9-]{3,})"),
               purrr::partial(.reformatAb, ab = !!ab),
               purrr::partial(separateSubunits, ab = !!ab, new_col = "subunit"),
               purrr::partial(splitMerge, ex = !!split_merge_str,
                              f = formatTCR, tcr = "greek_letter"),
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
    .warnIfColExists(df, new_col)

    df <- df %>%
        dplyr::mutate(!!new_col :=
                          do.call(paste, c(!!syms(id_cols), sep = sep))) %>%
        dplyr::group_by(!!sym(new_col)) %>%
        dplyr::mutate(n = n(),
                      i = row_number(),
                      !!new_col := ifelse(n == 1, !!sym(new_col),
                                     paste(!!sym(new_col), i, sep = sep)))

    if (isTRUE(warn) & max(df[, "n"]) > 1){
        warning("ID columns do not uniquely identify rows, row numbers added.")
    }

    df <- df %>% dplyr::select(-n, -i)

    return(df)
}


# dplyr version ----
#addID2 <- function(df, id_cols = c(Antigen, Study), new_col = ID,
#                  warn = TRUE){
#
#    #.warnIfColExists(df, new_col)
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
#' have new value NA
#'
#'@param ab (character(n)) A vector of strings to transform
#'@export
upperSquish <- function(ab){
    x <- toupper(gsub("(^[A-z0-9]+)[-\\. ]?([A-z0-9]+)$", "\\1\\2", ab))
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
                      upper_no_dash = upperSquish(greek_letter),
                      lower_no_dash = lowerNoDash(greek_letter),
                      dash_not_dot = dashNotDot(greek_letter))
    return(result)
}
