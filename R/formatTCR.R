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
#'
#'@author Helen Lindsay
#'@returns df, with an additional column, by default named "TCR_long",
#'giving the TCR name in long format to allow exact matching in the
#'gene_aliases table.
#'@examples
#'df <- data.frame(Antigen = c("TCRab", "TCR-ab", "TCRa/b"))
#'formatTCR(df, tcr = "Antigen")
formatTCR <- function(df, tcr="TCR", new_col="TCR_long"){
    .stopIfColExists(df, new_col)

    # Operate on the column to avoid worrying about temporary column names
    tcr_t <- unique(formatTCRv(df[[tcr]], new_col))
    colnames(tcr_t) <- c(tcr, new_col)

    # Merge result into original data.frame
    return(dplyr::full_join(df, tcr_t, by = tcr, multiple = "all"))

}


# formatTCRv ----
#'
#' Format TCR (vector)
#'
#' Used by formatTCR, accepts a vector (one column from a data.frame),
#' returns a long data.frame with a new column containing the TCR description
#'
#'@param tcr A vector of T cell receptor antibody names
#'@param new_col The name of the new column containing the TCR description
#'@importFrom dplyr all_of any_of
#'@keywords internal
#'@returns A data.frame, with column "TCR" being the original values of
#' parameter "tcr", and a second column "TCR_long" being the reformatted names
#' of each matching TCR subunits.
formatTCRv <- function(tcr, new_col="TCR_long"){

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
        tidyr::separate(.data$end, c("end", "nm"),
                        sep = "(?=[0-9])", fill = "right", extra = "merge") %>%
        # Create a column with the variable and joining prefixes
        tidyr::separate(.data$end, c("vj", "abgd"),
                        sep = "(?=[abgd])", extra = "merge") %>%
        dplyr::mutate(vj = stringr::str_replace_all(.data$vj,
                                                    var_join, names(var_join)),
                      vj = ifelse(.data$vj == "", "locus", .data$vj),
                      abgd = stringr::str_replace_all(.data$abgd,
                                                      l2w, names(l2w)),
                      !!new_col := gsub(" NA$", "",
                          paste("T cell receptor", .data$abgd,
                                .data$vj, .data$nm))) %>%
        dplyr::select(dplyr::all_of(c("TCR", new_col)))

    return(res)
}
