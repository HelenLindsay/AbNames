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
                                       list(!!dplyr::sym(ab)),
                                       !!dplyr::sym(temp_col)))
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
