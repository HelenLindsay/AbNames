# formatIg ----
#
#'Format immunoglobin antibody names
#'
#'@description
#'Adds a column to a data.frame containing immunoglobin names re-formatted to
#'match HGNC symbols.
#'
#'Assumes that Greek symbols have been converted to single letters.  Diversity
#'and joining regions are not implemented.  The final a/b for e.g. IgG2a is
#'dropped, as the specific modification isn't in protein ontology or HGNC.
#'The lambda light chain will only be matched in HGNC if a gene number is
#'provided.
#'
#'@param df A data.frame or tibble
#'@param ig (character(1), default "Antigen") Name of column containing antibody
#'names
#'@param new_col (character(1), default "Ig") Name of column containing
#'re-formatted immunoglobin names to add to df
#'@export
formatIg <- function(df, ig = "Antigen", new_col = "Ig"){
    .warnIfColExists(df, new_col)

    df <- dplyr::mutate(df, !!new_col := .formatIg(df[[ig]]))
    return(df)
}


#.formatIg ----
# Kappa light chain constant regions are IGKC
# There are several lambda constant regions - IGLC1, IGLC2 (predicted?)
# Heavy chain regions are named e.g. IGHM for IgM
.formatIg <- function(ig){
    # Note the IgG2 subclasses are covalent modifications, so I do not keep
    # the final A/B
    heavy <- .gsubNA("^Ig([ADEGM][0-4]?)[ABab]?[\\. ]?(Fc)?$",
                     sprintf("IGH%s", "\\1"), ig)
    kappa <- .gsubNA("^Ig .*k$", "IGKC", ig)
    lambda <- .gsubNA("^Ig.*l.*[0-9].*$", sprintf("IGLC%s", "\\1"), ig)
    result <- dplyr::coalesce(heavy, kappa, lambda)
    return(result)
}
