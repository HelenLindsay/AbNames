# searchTotalseq ----
#

#'Annotate gene IDs using totalseq data set
#'
#' Search for exact matches between antigen names in df and totalseq cocktails
#' or between catalogue numbers. Relies on column names matching totalseq
#' As the goal is to add gene IDs, clone is not used.  This function works by
#' successive left_joins and will add rows if columns HGNC_SYMBOL, HGNC_ID or
#' ENSEMBL_ID are already in x.
#'
#'@param x data.frame, containing at minimum columns named "Antigen"
#'@param cols list of columns to match sequentially in totalseq.  If not
#'specified, defaults to catalogue number then antigen then clone.
#'@export
searchTotalseq <- function(x, cols = NULL){
    if (is.null(cols)) {
        cols <- intersect(c("Cat_Number", "Antigen", "Clone"), colnames(x))
        cols <- as.list(cols)
    }


    id_cols <- c("ENSEMBL_ID", "HGNC_SYMBOL", "HGNC_ID")

    # If ID columns already exist, separate into filled and unfilled
    y <- head(x, 0)
    y <- x %>%
            dplyr::filter(if_any(dplyr::any_of(id_cols), ~ ! is.na(.x)))
    x <- dplyr::anti_join(x, y) %>%
        dplyr::select(-any_of(id_cols))

    utils::data("totalseq", envir = environment())

    # Select the gene information from totalseq
    ts <- totalseq %>%
        dplyr::select(Antigen, dplyr::all_of(c(unlist(cols), id_cols)))

    result <- left_join_any(x, ts, cols)
    return(dplyr::bind_rows(y, result))
}
