# searchTotalseq ----
#
# Search for exact matches between antigen names in df and totalseq cocktails
#
searchTotalseq <- function(df, ab = "Antigen"){
    # Should clone be used?
    utils::data("totalseq", envir = environment())

    # Select the gene information from totalseq
    ts <- totalseq %>%
        dplyr::select(Antigen, Ensembl_ID, Gene_Symbol)

    result <- dplyr::left_join(df, ts, by = structure("Antigen", names = ab))
    return(result)
}
