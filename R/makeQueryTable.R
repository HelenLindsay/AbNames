# makeQueryTable ----
#'@title makeQueryTable
#'@description Given a data.frame containing Antigen/Antibody names, reformat
#'the names into possible gene names and return the data.frame in long format.
#'@param df A data.frame or tibble
#'@param ab (character(1), default "Antigen") Name of the column in df
#' containing Antigen/Antibody names
#'@param id_cols (character(n)) Names of columns to paste to form ID column
#'@param funs A vector of string formatting functions to
#'@export
makeQueryTable <- function(df, ab = "Antigen",
                           id_cols = c("Antigen", "Study"), funs){

}



# addID ----
#'@title addID
#'@description Adds an ID column to a data frame
#'@param df A data.frame or tibble
#'@param id_cols (character(n)) Names of columns to paste to form ID column
#'@param new_col (character(1), default: "ID") Name of new ID column
#'@param warn (TRUE/FALSE, default: TRUE) If TRUE, warn if IDs are not unique
#'@return
#'@importFrom dplyr mutate group_by pull sym
addID <- function(df, id_cols = c("Antigen", "Study"), new_col = "ID",
                  warn = TRUE){

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


# removePrefix ----
#'@title removePrefix
#'@description Create a new column by removing a prefix from an existing
#'column.  The default prefix removes either "anti-" or "Anti-".  If prefix
#'is not found, value will be NA.
#'@param df A data.frame or tibble
#'@param ab (character(1), default "Antigen) Name of the column to remove
#'prefixes from
#'@param prefix (character(1)) A regular expression for matching prefixes
#'@param new_col (character(1), default NA Name of the column to add to df.
#'If NA, column ab is modified
removePrefix <- function(df, ab = "Antigen", prefix = "[Aa]nti-", new_col = NA){
    if (is.na(new_col)) new_col <- ab
    return(dplyr::mutate(df, !!new_col :=  gsub(prefix, "", !!sym(ab))))
}


# splitUnnest ----
#'@title splitUnnest
#'@description Convenience function for splitting a column at a delimiter,
#'unnesting (one row per value after splitting) and removing unnecessary
#'whitespace
#'@param df A data.frame or tibble
#'@param ab (character(1), default "Antigen) Name of the column to remove
#'prefixes from
#'@param pattern (character(1), default "[\\(\\)]") A regular expression for
#'using with strsplit. The default expression splits at "(" or ")".
#'@param new_col (character(1), default NA Name of the column to add to df.
#'If NA, column ab is modified
#'@importFrom tidyr unnest
#'@importFrom stringr str_squish
splitUnnest <- function(df, ab, pattern = "[\\(\\)]", new_col = NA){
    if (is.na(new_col)) new_col <- ab

    return(mutate(df, !!new_col :=
                      strsplit(!!sym(ab), pattern, perl = TRUE)) %>%
               unnest(cols = new_col) %>%
               mutate(!!new_col = str_squish(new_col)))
}


# separateComma ----
separateComma <- function(df, ab, delim = ", ", new_col = NA){
    if (is.na(new_col)) new_col <- ab

    return(mutate(df, !!new_col := strsplit(!!sym(ab), delim)) %>%
               unnest(cols = !!sym(new_col)))
}


# separateSubunits ----
separateSubunits <- function(df, ab){

}


# replaceGreek ----
replaceGreek <- function(df, ab){

}

