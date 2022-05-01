# formatTCR ----
#
#'Convert T-cell receptor antigen names
#'
#'Convert T-cell antigen names to long format used in HGNC gene descriptions,
#'e.g. "T cell receptor gamma variable 24.  Greek symbols and words should be
#'converted to letters before using this function
#'
#'@param tcr (character(n)) List of unformatted TCR names
#'@importFrom stringr str_replace_all
formatTCR <- function(tcr){

    letter2word <- structure(c("alpha ", "beta ", "gamma ", "delta "),
                             names = c("a", "b", "c", "d"))

    df <- data.frame(tcr = tcr,
                     tcr_end = gsub("TCR\\s?", "", gsub("\\.", "-", tcr))) %>%
        splitUnnest("tcr_end", split = "\\/") %>%
        splitUnnest("tcr_end", split = "-(?=J)")





    variable <- .gsubNA("^.*V([abgd][0-9-]+).*$", "\\1", tcr_end)
    variable <- gsub("-$", "", variable_n)

}
