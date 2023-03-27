# matchToCiteseq -----
# Probably want to standardise name via aliases table.....
# At the moment relies on column names matching
# Rename programmatically?
# Match to the curated citeseq data.set
# nm could be either data.frame or single vector of names
# Eventually, matching just by names may be enough
# To do: check types of matching columns?
# To do: add a column indicating group size?
#
#' Standardise names using the citeseq data set
#'
#'@description Groups antibody names and selects the most frequent name.  By
#'default, the citeseq dataset is grouped by Antigen, Clone,
#'Cat_Number (catalog number) and ALT_ID.
#'@param x A data.frame or tibble containing a column "Antigen" to match to the
#' citeseq data set
#'@param cols (character(n), default NULL) Optional additional columns to use
#'for matching in citeseq data set.  ID columns from the citeseq dataset are
#'used by default.  Antibodies are grouped by a match in any of the ID columns.
#'@param verbose Report which columns are used for matching? (default: TRUE)
#'@param check.matches Set to TRUE to check which entries were matched
#'(default FALSE)
#'@param ... Options passed on to getCommonName
#'@returns If dry.run is FALSE, the data.frame x is returned with an additional
#'two columns "Antigen_std", containing the standardised name, and "n_matched",
#'the number of entries used for deriving the standard name.  The minimum
#'value of "n_matched" is 1 and means that no matches were found.
#'
#'If dry.run is TRUE, a grouped tibble is returned, containing rows from x
#'grouped with rows from the citeseq data set that were matched.  This can be
#'used for checking whether matches are correct and identifying what was matched
#'incorrectly if they are incorrect.
#'@export
matchToCiteseq <- function(x, cols = NULL, verbose = TRUE,
                           check.matches = FALSE, ...){
    # Load citeseq data set
    data_env <- new.env(parent = emptyenv())
    utils::data("citeseq", envir = data_env, package = "AbNames")
    citeseq <- data_env[["citeseq"]]

    # Set up names of columns for matching
    keep_cols <- intersect(cols, colnames(citeseq))

    # Check that inputs are okay
    .check_matchToCiteseq_inputs(x, cols, keep_cols)

    cn_x <- colnames(x)

    # If keep_cols are specified, make sure Antigen is included
    if (! is.null(keep_cols)){
        keep_cols <- unique(c(keep_cols, "Antigen"))
        msg_cols <- toString(keep_cols)
    } else {
        # If columns are not specified, check which of the default columns
        # in getCommonName are present in x
        default_cols <- c("Antigen", "Cat_Number", "Clone", "ALT_ID")
        msg_cols <- toString(intersect(colnames(x), default_cols))
    }

    # Report which columns are used to match new data to cite seq data
    msg <- "Matching new data to citeseq data using columns:\n%s"
    message(sprintf(msg, msg_cols))

    # Add temporary ID column and add new data to citeseq data set
    id <- .tempColName(x, nm = "ID")
    x <- x %>% dplyr::mutate(!!id := "KEEPME")
    x <- dplyr::bind_rows(x, citeseq)

    # If keep_cols is NULL, the default column set is used for matching
    x <- getCommonName(x, cols = keep_cols, ab = "Antigen",
                       fill_col = "Antigen_std", keep = TRUE, ...)

    #### HERE
    if (isTRUE(check.matches)){
        stop("IMPLEMENT check.matches")
    }

    # Remove ID column, select original columns
    x <- x %>% dplyr::filter(!!sym(id) == "KEEPME") %>%
        dplyr::select(all_of(c("Antigen_std", cn_x, "n_matched")))

    return(x)
}


# .matchToCiteseq_input_checks ----
.check_matchToCiteseq_inputs <- function(x, cols, keep_cols){
    if (! "data.frame" %in% class(x)) {
        stop("x should be a data.frame or tibble")
    }

    if (! "Antigen" %in% colnames(x)){
        stop("x must contain a column 'Antigen'")
    }

    if (! identical(keep_cols, cols)){
        warning("cols should match columns in cite seq data set")
    }
}


# Data for using in tests?
## Some non-standard names, including duplicates
#adt_nms <- c("KLRG1 (MAFA)", "KLRG1", "CD3 (CD3E)", "HLA.A.B.C",
#             "HLA-A/B/C", "NKAT2", "CD66a.c.e", "CD66a_c_e",
#             "CD11a/CD18 (LFA-1)")     # NKAT2 = CD158b

