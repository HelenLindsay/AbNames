# searchTotalseq ----
#

#'Annotate gene IDs using totalseq data set
#'
#' Search for exact matches between antigen names in df and totalseq cocktails
#' or between catalogue numbers. Relies on column names matching totalseq
#' As the goal is to add gene IDs, clone is not used.
#'
#'@param x data.frame, containing at minimum columns named "Antigen"
#'@param cols list of columns to match sequentially in totalseq.  If not
#'specified, defaults to catalogue number then antigen then clone.
#'@export
searchTotalseq <- function(x, cols = NULL){
    if (is.null(cols)) {
        cols <- list("Cat_Number", "Antigen", "Clone")
    }

    # Select the gene information from totalseq
    ts <- totalseq %>%
        dplyr::select(Antigen, ENSEMBL_ID, HGNC_SYMBOL, HGNC_ID,
                      any_of(unlist(cols)))

    utils::data("totalseq", envir = environment())
    return(left_join_any(x, ts, cols))
}
