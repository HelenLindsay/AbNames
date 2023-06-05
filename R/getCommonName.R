# getCommonName ----
#'Find most common antibody name
#'
#'Find most common name by matching any/all of Antigen name, clone, gene
#'identifiers
#'
#'Be careful of exceptions, e.g. Cat_Number == "custom made"
#'
#' @param x data.frame for finding common name
#' @param cols Columns for grouping, default: NULL means that the columns
#' Antigen, Cat_Number, Clone, ALT_ID are used.  Note that gene identifiers
#' are not used by default as some antibodies target particular isoforms.
#' @param ab Column containing antibody names for standardising,
#'default "Antigen"
#' @param fill_col column to add
#' @param n_matched name of column to add with number of matching results,
#' default "n_matched"
#' @param ignore A named list, where names match colnames of x and values are
#' regular expressions matching values in x that should not be used for
#' grouping. For example, entries where the catalogue number
#' is "Custom made" should not be used for matching.  Alternatively, set to TRUE
#' to use the default ignore patterns.
#' @param id Name of column containing ID, will be "KEEPME" for the query data
#' if called via matchToCiteseq.  If present, citeseq data will be subset after
#' matching to include only one row per Antigen/Study combination before
#' standardising names.
#' @param ... pass keep = TRUE for keeping grouping columns for debugging
#' @param verbose Print information about which columns are used for grouping?
#' Default: TRUE
#' @author Helen Lindsay
#' @returns x, with additional columns fill_col (default name "Antigen_std")
#' with the standardised Antigen name; "group" containing the grouping indices;
#' and n_matched (default name "n_matched"), containing the number of rows
#' that were used to generate the standard name.
#' @examples
#' df <- data.frame(Antigen = c("B220", "B220 (CD45R)",
#'                              "C-kit (CD117)", "C-KIT"),
#'                  HGNC_ID = rep(c("HGNC:9666", "HGNC:6342"), each = 2))
#'
#'# Here we match Antigens using HGNC_ID.  Note that we have to explicitly
#'# specify the matching columns.  By default, HGNC_ID is not used for matching
#'# to avoid matching isoforms.
#'getCommonName(df, cols = "HGNC_ID")
getCommonName <- function(x, cols=NULL, ab="Antigen",
                          fill_col="Antigen_std", n_matched="n_matched",
                          ignore=TRUE, id=NULL, verbose=TRUE, ...){

    .check_getCommonName(colnames(x), ab, fill_col, n_matched, cols)

    if (isTRUE(ignore)){
        ignore <- list(Cat_Number = "[Cc]ustom", Clone = CLONE_DUPS()$Clone)
    }

    if (is.null(cols)) {
        cols <- intersect(c("Antigen", "Cat_Number", "Clone", "ALT_ID"),
                          colnames(x))
    }

    # Initialise fill col to Antigen up until first bracket
    # (brackets are usually alternative names)
    x <- dplyr::mutate(x, !! fill_col := gsub(" ?\\(.*", "", !!sym(ab)))

    tmp_grp <- .tempColName(x, nm="group") # Don't overwrite existing "group"

    if (isTRUE(verbose)){
        message(sprintf("Grouping data using columns:\n%s", toString(cols)))
    }

    # Group by any e.g. catalogue number or exact match to antigen
    x <- group_by_any(x, groups=cols, new_col=tmp_grp, ignore=ignore)

    # Make sure that one Antigen is only counted once per study
    if (! is.null(ID)){
        query <- x %>% dplyr::filter(.data[[id]] == "KEEPME")
        x <- x %>%
            dplyr::filter(.data[[id]] != "KEEPME") %>%
            dplyr::group_by(.data[[tmp_grp]], .data$Study) %>%
            dplyr::slice(1) %>%
            dplyr::bind_rows(query) %>%
            dplyr::group_by(.data[[tmp_grp]])
    }

    x <- x %>% dplyr::mutate(!!n_matched := dplyr::n())

    # Fill with most common value
    x <- fillByGroup(x, group=tmp_grp, method="all",
                     multiple="mode", fill=fill_col)

    return(x)
}


# .check_getCommonName input checks ----
.check_getCommonName <- function(col_nms, ab, fill_col, n_matched, cols){
    if (! ab %in% col_nms){
        stop("Parameter 'ab' must be a column name in x (default 'Antigen')")
    }

    msg <- "Column %s already exists in data.frame"
    if (fill_col %in% col_nms){ stop(sprintf(msg, fill_col)) }

    if (n_matched %in% col_nms){ stop(sprintf(msg, n_matched)) }

    if (! all(cols %in% col_nms)){
        stop("Parameter 'cols' must be a vector of column names in x,",
             "which will be used for matching Antigen names")
    }
}
