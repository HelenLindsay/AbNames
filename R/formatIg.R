#
# Kappa light chain constant regions are IGKC
# There are several lambda constant regions - IGLC1, IGLC2 (predicted?)
# Heavy chain regions are named e.g. IGHM for IgM
#
# Assume that Greek symbols have been converted to single letters.  I have not
# implemented diversity or joining regions as I haven't seen CITE-seq antibodies
# specifically for these.  Final AB for IgG2 is dropped, the specific
# modification isn't in protein ontology or HGNC
formatIg <- function(df, ig = "Antigen", new_col = "Ig"){
    warnIfColExists(df, new_col)

    df <- dplyr::mutate(df, !!new_col := .formatIg(df[[ig]]))
    return(df)
}


.formatIg <- function(ig){
    # Note the IgG2 subclasses are covalent modifications, so I do not keep
    # the final A/B
    heavy <- .gsubNA("^Ig([ADEGM][0-4]?)[ABab]?[\\. ]?(Fc)?$",
                     sprintf("IGH%s", "\\1"), ig)
    kappa <- .gsubNA("^Ig .*k$", "IGKC", ig)
    lambda <- .gsubNA("^Ig .*l$", "IGLC@", ig) # Match the whole cluster
    result <- dplyr::coalesce(heavy, kappa, lambda)
    return(result)
}
