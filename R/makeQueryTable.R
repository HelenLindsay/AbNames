# makeQueryTable ----
#'Separate Antibody/Antigen names into component parts
#'
#'Given a data.frame containing Antigen/Antibody names, reformat
#'the names into possible gene names and return the data.frame in long format.
#'
#'@param df A data.frame or tibble
#'@param funs A vector of string formatting functions to apply to df.  Each
#'should take a data.frame as a first argument and return a data.frame or
#'tibble, and optional arguments should be filled in e.g. with purrr::partial.
#'@param control_col (character(1), default: NA)  Optional name of a logical
#'column indicating whether an antibody is an isotype control.  If present,
#'controls will be removed to avoid spurious matches as they usually do not
#'react against human genes.
#'@param query_cols (character(n), default: NA) IF NA ASSUME THAT ALL COLUMNS
#'ADDED ARE QUERY_COLS?
#'@importFrom tidyr pivot_longer
#'@export
makeQueryTable <- function(df, funs, query_cols = NA, control_col = NA){
    # Remove controls
    if (! is.na(control_col)) { df <- dplyr::filter(df, ! (!!sym(control_cl))) }

    cn <- colnames(df)

    df <- magrittr::freduce(df, funs)

    if (is.na(query_cols)) query_cols <- setdiff(colnames(df), cn)

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
defaultQuery <- function(ab = "Antigen", id_cols = c("Antigen", "Study")){

    # Default transformation sequence for making query table ----

    # List of functions, column to act on and column to create?
    # Check that all columns to act on are either in the table or created
    # earlier in the pipeline
    # function, ab, new_col, args?
    query_cols <- c(ab)

    qry = list(purrr::partial(addID, id_cols = id_cols),
               purrr::partial(gsubAb, ab = ab), # Remove A/antis
               purrr::partial(gsubAb, ab = ab, pattern = "\\sRecombinant"),
               purrr::partial(splitUnnest, ab = ab), # Brackets
               purrr::partial(splitUnnest, ab = ab, split = ", "),  # Commas
               # / _ or . if at least 3 on each side and not TCR
               purrr::partial(splitUnnest, exclude = "TCR",
                        split = "(?<=[A-z0-9-]{3})[\\/_\\.](?=[A-z0-9-]{3,})")
    )

    return(list(query_funs = qry, query_cols = query_cols))
}


# addID ----
#'add ID column
#'
#'Adds an ID column to a data frame
#'
#'@param df A data.frame or tibble
#'@param id_cols (character(n)) Names of columns to paste to form ID column
#'@param new_col (character(1), default: "ID") Name of new ID column
#'@param warn (TRUE/FALSE, default: TRUE) If TRUE, warn if IDs are not unique
#'
#'@return df with an extra ID column
#'
#'@importFrom dplyr mutate group_by pull all_of n syms
addID <- function(df, id_cols = c("Antigen", "Study"), new_col = "ID",
                  warn = TRUE){

    # TO DO:
    # should - be replaced by _ in id columns to make unsplitting easier?

    # Check id_cols exist in df and new_col does not
    if (! all(id_cols %in% colnames(df))){ stop("All id_cols must be in df") }
    .warnIfColExists(df, new_col)

    df <- df %>%
        dplyr::mutate(!!new_col :=
                          do.call(paste, c(!!syms(id_cols), sep = "-")))

    if (isTRUE(warn)){

        # Check that the ID column uniquely identifies rows
        max_group_n <- group_by(df, !!sym(new_col)) %>%
            mutate(n = n()) %>%
            pull(n) %>%
            max()

        if (max_group_n > 1){
            warning("There are duplicate values in ID column")
        }
    }

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
    return(.noDups(x, ab))
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


# splitUnnest ----
#'Split a column and create one row per entry
#'
#'Convenience function for splitting a column at a delimiter,
#'unnesting (one row per value after splitting) and removing unnecessary
#'whitespace.  Default is to split at brackets.  Returns a tibble tbl_df.
#'
#'@param df A data.frame or tibble
#'@param ab (character(1), default "Antigen) Name of the column to remove
#'prefixes from
#'@param split (character(1), default "[\\(\\)]") A regular expression for
#'using with strsplit. The default expression splits at "(" or ")".
#'@param new_col (character(1), default NA Name of the column to add to df.
#'If NA, column ab is modified.
#'@param exclude (character(1), default NA) a regex - do not split if ab
#'matches.
#'
#'@examples
#'df <- data.frame(Antigen = c("CD279 (PD-1)", "Mac-2 (Galectin-3)"))
#'splitUnnest(df, "Antigen", new_col = "Split")
#'
#'@importFrom tidyr unnest
#'@importFrom stringr str_squish
#'@importFrom dplyr sym select rename all_of
#'@export
splitUnnest <- function(df, ab = "Antigen", split = "[\\(\\)]", new_col = NA,
                        exclude = NA){

    temp_col <- .tempColName(df)

    df <- dplyr::mutate(df, !!temp_col :=
                            strsplit(!!dplyr::sym(ab), split, perl = TRUE))

    if (! is.na(exclude)){
        df <- dplyr::mutate(df, !!temp_col :=
                ifelse(grepl(exclude, !!dplyr::sym(ab)),
                       !!dplyr::sym(ab), !!dplyr::sym(temp_col)))
    }

    df <- tidyr::unnest(df, cols = all_of(temp_col)) %>%
        dplyr::mutate(!!temp_col := str_squish(!!sym(temp_col)))

    # If temp_col should overwrite ab, remove ab and rename new_col
    if (is.na(new_col)){
        df <- df %>% dplyr::select(-!!dplyr::sym(ab))
        new_col <- ab
    }

    df <- df %>% dplyr::rename(!!new_col := !!dplyr::sym(temp_col))
    return(df)
}

