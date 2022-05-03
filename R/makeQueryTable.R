# makeQueryTable ----
#'Separate Antibody/Antigen names into component parts
#'
#'Given a data.frame containing Antigen/Antibody names, reformat
#'the names into possible gene names and return the data.frame in long format.
#'
#'@param df A data.frame or tibble
#'@param ab (character(1), default "Antigen") Name of the column in df
#' containing Antigen/Antibody names
#'@param id_cols (character(n)) Names of columns to paste to form ID column
#'@param funs A vector of string formatting functions to
#'
#'@export
makeQueryTable <- function(df, ab = "Antigen",
                           id_cols = c("Antigen", "Study"), funs){

}


defaultQry <- function(){

    # Default transformation sequence for making query table ----

    # List of functions, column to act on and column to create?
    # Check that all columns to act on are either in the table or created
    # earlier in the pipeline
    # function, ab, new_col, args?
    stopifnot(requireNamespace("pryr"))

    qry = list(addID,
               gsubAb(), # Remove A/anti
               pryr::partial(gsubAb, pattern = "\\sRecombinant"),
               splitUnnest(), # Brackets
               splitUnnest(split = ", "),  # Commas
               # / _ or . if at least 3 on each side and not TCR
               pryr::partial(splitUnnest, exclude = "TCR",
                        split = "(?<=[A-z0-9-]{3})[\\/_\\.](?=[A-z0-9-]{3,})")

    )
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
#'@return
#'
#'@importFrom dplyr mutate group_by pull all_of n syms
addID <- function(df, id_cols = c("Antigen", "Study"), new_col = "ID",
                  warn = TRUE){

    # Check id_cols exist in df and new_col does not
    if (! all(id_cols %in% colnames(df))){ stop("All id_cols must be in df") }
    .warnIfColExists(df, new_col)

    df <- mutate(df, !!new_col := do.call(paste, c(!!syms(id_cols), sep = "-")))

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
#'@param new_col (character(1), default NA Name of the column to add to df.
#'If NA, column ab is modified
#'
gsubAb <- function(df, ab = "Antigen", pattern = "[Aa]nti-", replacement = "",
                   exclude = NA, restrict = NA, new_col = NA){
    if (is.na(new_col)) new_col <- ab
    df <- dplyr::mutate(df, !!new_col :=  gsub(pattern, replacement, !!sym(ab)))

    # Restrict would have to be a filter expression, e.g. a particular study
    # Would need to use a temp column as in splitUnnest
    # Is this actually important?

    #@param exclude (character(1), default NA) a regex - do not split if ab
    #matches.
    #@param restrict (character(1), default NA) a regex replace only if restrict
    # matches in column ab

    return(df)
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




