# searchTotalseq ----
#
#' Annotate antibodies with gene IDs using totalseq data
#'
#' Search for exact matches between antigen names, clones, or catalogue numbers
#' in data.frame x and totalseq cocktails data set.  Relies on column names in
#' matching column names in totalseq (Antigen, Cat_Number and/or Clone).  In our
#' experience, the antibody clone almost always uniquely identifies an antibody,
#' but there is at least one exception so we recommend caution and manual
#' checking if matching by only Clone.
#'
#' This function works by repeated left_joins and can add rows.  Only rows
#' that do not already have identifiers are joined with the totalseq data.  This
#' function does not check if existing identifiers match totalseq identifiers.
#' If antigens are matched by Cat_Number or Clone and Antigen is NA, Antigen
#' will be filled.
#'
#'@param x data.frame, containing at minimum columns named "Antigen"
#'@param cols list of columns to match in totalseq data.  If not
#'specified, defaults to catalogue number, antigen and clone.
#'@author Helen Lindsay
#'@export
#'@returns Data.frame x with additional identifier columns from joining x
#' with the totalseq data set.
#'@examples
#'searchTotalseq(data.frame(Cat_Number = c("306623", "124235", "331609")))
searchTotalseq <- function(x, cols=NULL){
    # Which columns should be used for matching?
    if (is.null(cols)) {
        cols <- intersect(c("Cat_Number", "Antigen", "Clone"), colnames(x))
    }

    # These are the names of ID columns contained in totalseq data
    id_cols <- c("ENSEMBL_ID", "HGNC_SYMBOL", "HGNC_ID", "ALT_ID")

    # If ID columns already exist, separate into filled and unfilled
    y <- x %>%
            dplyr::filter(if_any(dplyr::any_of(id_cols), ~ ! is.na(.x)))

    if (! any(id_cols %in% colnames(x))) { y <- utils::head(x, 0) }

    x <- dplyr::anti_join(x, y, by=colnames(x)) %>%
        dplyr::select(-dplyr::any_of(id_cols))

    # Load totalseq data set
    data_env <- new.env(parent=emptyenv())
    utils::data("totalseq", envir=data_env, package="AbNames")
    totalseq <- data_env[["totalseq"]]

    # Select the gene information from totalseq
    ts <- totalseq %>%
        dplyr::select(dplyr::all_of(c("Antigen", cols, id_cols)))

    result <- left_join_any(x, ts, cols)
    result <- dplyr::bind_rows(result, y)

    # If Antigen is missing but Cat_Number or Clone is matched, patch Antigen
    fill_cols <- intersect(cols, c("Cat_Number", "Clone"))
    if (any(is.na(result$Antigen)) & length(fill_cols) >= 1){
        ts <- ts %>%
            filter_by_union(result %>% dplyr::filter(is.na(.data$Antigen))) %>%
            unique()
        result <- dplyr::rows_patch(result, ts, unmatched = "ignore",
                                    by=fill_cols)
    }

    # For entries matched by Cat_Number or Clone, fill any matching Antigens
    result <- result %>%
        dplyr::group_by(Antigen) %>%
        tidyr::fill(dplyr::any_of(id_cols), .direction="updown") %>%
        dplyr::ungroup()

    return(result)
}
