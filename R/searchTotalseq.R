# searchTotalseq ----
#
#'Annotate gene IDs using totalseq data set
#'
#' Search for exact matches between antigen names or catalogue numbers in
#' data.frame x and totalseq cocktails data set.  Relies on column names in
#' matching column names in totalseq.
#' This function works by repeated left_joins and can add rows.
#'
#'@param x data.frame, containing at minimum columns named "Antigen"
#'@param cols list of columns to match in TotalSeq.  If not
#'specified, defaults to catalogue number, antigen and clone.
#'@export
searchTotalseq <- function(x, cols = NULL){
    # Which columns should be used for matching?
    if (is.null(cols)) {
        cols <- intersect(c("Cat_Number", "Antigen", "Clone"), colnames(x))
        cols <- as.list(cols)
    }

    # These are the names of ID columns contained in totalseq data
    id_cols <- c("ENSEMBL_ID", "HGNC_SYMBOL", "HGNC_ID", "ALT_ID")

    # If ID columns already exist, separate into filled and unfilled
    y <- utils::head(x, 0)
    y <- x %>%
            dplyr::filter(if_any(dplyr::any_of(id_cols), ~ ! is.na(.x)))
    x <- dplyr::anti_join(x, y) %>%
        dplyr::select(-dplyr::any_of(id_cols))

    utils::data("totalseq", envir = environment())

    # Select the gene information from totalseq
    ts <- totalseq %>%
        dplyr::select(dplyr::all_of(c("Antigen", unlist(cols), id_cols)))

    result <- left_join_any(x, ts, cols)
    return(dplyr::bind_rows(y, result))
}
