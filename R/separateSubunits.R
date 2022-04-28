# separateSubunits ----
#' Separate multi-subunit protein names
#'
#'@description
#' Separate names of antibodies against multi-subunit proteins
#' e.g. CD235ab, CD66ace into one subunit per row.
#'
#'Two subunit patterns are considered.  For the first, subunits are lower case
#'letters and the gene name has no separator, e.g. CD66ace is composed of
#'subunits CD66a, CD66b and CD66c.  For the second pattern, subunits are written
#'with uppercase letters and are separated with a "-", e.g. HLA-A/C/E is
#'composed of subunits HLA-A, HLA-C and HLA-E.  Both patterns require at least
#'at least 2 capital letters or numbers followed by at least 2 possible
#'subunits. There may be a separator between the groups and/or between
#'the lower case letters.  At present, the between group separators
#'are -, . and space, and the between subunit separators are / and .
#'
#'Subunits should be converted from Greek symbols before applying this function.
#'
#'
#'@param df A data.frame or tibble
#'@param ab (character(1), default "Antigen) Name of the column containing
#'antibody names
#'@param new_col (default: subunit) Name of new column containing guesses for
#'single subunit names
#'@importFrom dplyr case_when
separateSubunits <- function(df, ab, new_col = "subunit"){
    tmp_cn <- .tempColName(df, n = 4)

    df <- df %>%
        dplyr::mutate(
            # Select the subunit ends:
            # Pattern 1: At least 2 capital letters/numbers, optional separator,
            # then at least 2 lowercase with optional separator
            tmp_cn[1] :=  gsub("^[A-Z0-9]{2,}[-\\. ]?([a-z\\/\\.]{2,})",
                               "\\1", !!sym(ab)),

            # Pattern 2: At least 2 capital letters/numbers, then -, then
            # at least 2 uppercase letters or numbers with optional / or ,
            tmp_cn[2] := gsub("^[A-Z0-9]{2,}-([A-Z0-9\\/,]{2,})",
                              "\\1", !!sym(ab)),

            # Set to NA if neither pattern was matched, otherwise pick the match
            tmp_cn[3] := case_when(!!sym(tmp_cn[1]) == ab &
                                        !!sym(tmp_cn[2]) == ab ~ NA_character_,
                                    !!sym(tmp_cn[1]) != ab ~ !!sym(tmp_cn[1]),
                                    !!sym(tmp_cn[2]) != ab ~ !!sym(tmp_cn[2]),
                                    TRUE ~ "FIXME")
        )

}


#.separateSubunits ----
#'
#'@param df
#'@param ab
#'@param pattern
#'@param join_pattern sprintf pattern for joining t1 (start) and t2 (end)
#'@param t1
#'@param t2
.separateSubunits <- function(df, ab, new_col, pattern, join_pattern, t1, t2){
    df <- df %>%
        dplyr::mutate(t2 := gsubNA(pattern, "\\1", !!sym(ab)),
                      t1 := stringr::str_replace(!!sym(ab), !!sym(t2), ""),
                      t1 := gsub("[- ]$", "", !!sym(t1)),
                      t2 := strsplit(gsub("[,/\\.]", "", !!sym(t2)), "")) %>%
        tidyr::unnest(t2) %>%
        dplyr::mutate(!!new_col :=
                          sprintf(join_pattern, !!sym(t1), !!sym(t2))) %>%
        dplyr::select(-!!sym(t1), -!!sym(t2))
}



