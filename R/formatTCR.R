# formatTCR ----
#
#'Convert T-cell receptor antigen names
#'
#'Convert T-cell antigen names to long format used in HGNC gene descriptions,
#'e.g. "T cell receptor gamma variable 24.  Greek symbols and words should be
#'converted to letters before using this function
#'
#'@param df data.frame or tibble
#'@param tcr (character(1)) Name of column in df containing TCR names
#'@param new_col (character(1)) Name of column to
#'@return df in long form with one subunit per row and a new column "new_col"
#'containing the TCR gene description query string
#'@importFrom stringr str_replace_all
#'@importFrom tidyr separate
#'@importFrom dplyr rename
#'@export
formatTCR <- function(df, tcr = "TCR", new_col = "TCR_long"){
    .warnIfColExists(df, new_col)

    # Operate on the column to avoid worrying about temporary column names
    tcr_t <- formatTCRv(df[[tcr]], new_col) %>%
        unique() %>%
        dplyr::rename(!!tcr := TCR)

    # Merge result into original data.frame
    return(dplyr::full_join(df, tcr_t))

}


# formatTCRv ----
#'
#' Format TCR (vector)
#'
#' Used by formatTCR, accepts a vector (one column from a data.frame),
#' returns a long data.frame with a new column containing the TCR description
formatTCRv <- function(tcr, new_col = "TCR_long"){

    # Substitute Greek letters for Greek words for cases:
    # e.g. a (single letter) or "Va7-2" (preceded by V/J, followed by numbers)
    l2w <- structure(c("alpha", "beta", "gamma", "delta"),
                     names = sprintf("^%s[0-9-]*$", c("a", "b", "g", "d")))

    var_join <- structure(c("variable", "joining", "constant"),
                          names = c("V", "J", "C"))

    # Substitute . or _ with dash, remove prefix TCR with optional space or -
    res <- data.frame(TCR = tcr,
                      end = gsub("TCR[\\s-]?", "",
                                 gsub("[\\._]", "-", tcr))) %>%
        splitUnnest("end", split = "\\/") %>%
        splitUnnest("end", split = "-(?=J)") %>% # - if followed by J
        splitUnnest("end", split = "(?<=a)(?=b)|(?<=g)(?=d)") %>% # ab or gd
        tidyr::separate(end, c("end", "nm"),
                        sep = "(?=[0-9])", fill = "right", extra = "merge") %>%
        # Create a column with the variable and joining prefixes
        tidyr::separate(end, c("vj", "abgd"),
                        sep = "(?=[abgd])", extra = "merge") %>%
        dplyr::mutate(vj = stringr::str_replace_all(vj,
                                                    var_join, names(var_join)),
                      vj = ifelse(vj == "", "locus", vj),
                      abgd = stringr::str_replace_all(abgd, l2w, names(l2w)),
                      !!new_col := gsub(" NA$", "",
                          paste("T cell receptor", abgd, vj, nm))) %>%
        dplyr::select(TCR, !!sym(new_col))

    return(res)
}
