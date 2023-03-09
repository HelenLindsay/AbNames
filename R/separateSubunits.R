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
#'At present user-supplied regex patterns are not supported
#'
#'@param df A data.frame or tibble
#'@param ab (character(1), default "Antigen) Name of the column containing
#'antibody names
#'@param new_col (default: subunit) Name of new column containing guesses for
#'single subunit names
#'@importFrom dplyr case_when coalesce
#'@export
separateSubunits <- function(df, ab = "Antigen", new_col = "subunit"){
    tmp <- .tempColName(df, n = 4)

    # Pattern 1: At least 2 capital letters/numbers, optional separator,
    # then at least 2 lowercase with optional separator, not more than 8
    p1 <- "^[A-Z0-9]{2,}[-\\. ]?([a-z\\/\\.]{2,6})$"

    # Pattern 2: At least 2 capital letters/numbers, then - or ., then
    # at least 2 uppercase letters or numbers with optional / . or ,
    p2 <- "^[A-Z0-9]{2,}[-\\.]([A-Z\\/,\\.]{2,6})$"

    df <- .separateSubunits(df, ab, tmp[3], p1, "%s%s", tmp[1], tmp[2])
    df <- .separateSubunits(df, ab, tmp[4], p2, "%s-%s", tmp[1], tmp[2])

    # Merge the two patterns
    df <- df %>%
        dplyr::mutate(!!new_col :=
                          dplyr::coalesce(!!sym(tmp[3]), !!sym(tmp[4]))) %>%
        dplyr::select(-!!sym(tmp[3]), -!!sym(tmp[4]))

    return(df)
}


#.separateSubunits ----
#'
#'Separate multi-subunit protein names into one name per row
#'
#'Separates a column at a provided pattern, separates the second second into
#'individual letters and joins segments with a provided joining pattern
#'
#'@param df A data.frame or tibble
#'@param ab The name of the column containing names to split
#'@param new_col The name of the new column containing split names
#'@param pattern The regex pattern to use for splitting
#'@param join_pattern sprintf pattern for joining t1 (start) and t2 (end)
#'@param t1 first temporary column name
#'@param t2 second temporary column name
#'@importFrom rlang :=
#'@importFrom dplyr all_of
.separateSubunits <- function(df, ab, new_col, pattern, join_pattern, t1, t2){
    #  If there are any duplicated characters, it's probably not a subunit
    no_dup <- function(x){
        n_char <- lengths(x)
        # CHECK - CHANGED FROM SAPPLY
        n_unq <- vapply(x, function(y) length(unique(y)), integer(1))
        return(n_char == n_unq)
    }

    df <- df %>%
        # Select the end pattern, remove end pattern and space to get start,
        # split end pattern (subunits).  If any subunits repeated once split,
        # set to NA because it's probably incorrect
        dplyr::mutate(!!t2 := .gsubNA(pattern, "\\1", !!sym(ab)),
                      !!t1 := stringr::str_replace(!!sym(ab), !!sym(t2), ""),
                      !!t1 := gsub("[- ]$", "", !!sym(t1)),
                      !!t2 := strsplit(gsub("[,/\\.]", "", !!sym(t2)), ""),
                      !!t2 := ifelse(no_dup(!!sym(t2)), !!sym(t2), list(NA))) %>%

        # Unnest to give one row per subunit
        tidyr::unnest(all_of(t2)) %>%
        # Join subunits
        dplyr::mutate(!!new_col :=
                          .printf(join_pattern, !!sym(t1), !!sym(t2))) %>%
        # Remove temporary columns
        dplyr::select(-!!sym(t1), -!!sym(t2))
    return(df)
}


# .checkSubunitMatches ----
#
#' Check if all subunits of a multi-subunit protein are matched
#
#' Internal AbNames function for removing spurious matches caused by guessing
#' subunit names.  Column names are hard-coded and expected to match the default
#' pipeline.  Returns a table containing gene name matches with spurious matches
#' removed.
#'
#'@param df A data frame containing matches of gene names in a database
#'@param query_df A data.frame containing potential gene/protein names, one
# per row
#'@importFrom dplyr anti_join
#'@importFrom rlang .data
.checkSubunitMatches <- function(df, query_df){
    nms <- c("subunit", "TCR_long")

    nsubunits <- query_df %>%
        dplyr::filter(.data$name %in% nms) %>%
        dplyr::group_by(.data$ID, .data$name) %>%
        dplyr::summarise(nexpected = dplyr::n())

    incomplete <- df %>%
        dplyr::filter(.data$ID %in% nsubunits$ID, .data$name %in% nms) %>%
        dplyr::group_by(.data$ID, .data$name) %>%
        dplyr::summarise(nmatched = dplyr::n()) %>%
        dplyr::full_join(nsubunits, by = c("ID", "name")) %>%
        dplyr::filter(! .data$nmatched == .data$nexpected) %>%
        dplyr::select(-.data$nmatched, -.data$nexpected)

    result <- df %>%
        dplyr::anti_join(incomplete, by = c("ID", "name"))

    return(result)
}

