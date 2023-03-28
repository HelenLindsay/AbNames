# sharedSubstr ----
#
#' Group based on shared substrings
#'
#' Takes a data.frame containing a column of words (e.g. Antibody aliases)
#' and a corresponding column of grouping ids
#' and returns a numeric vector indicating groupings where members of a
#' group share a word with at least one other member of the group.
#'
#' The intention is that if an antibody has alternative names and is sometimes
#' called by both, e.g. CD274 (B7-H1), we would like to group all entries
#' matching CD274 with all entries matching B7-H1. Intended to be called as
#' part of a querying pipeline. Names must already be broken into words, e.g.
#' for CD274 (B7-H1), CD274 and B7-H1 should appear in separate rows with the
#' same ID, e.g. by first calling [splitUnnest()].
#'
#' Although this function can be useful for matching antibody names, in our
#' experience manual checking of the results is required.
#'
#'@param df A query data.frame, e.g. created by makeQueryTable
#'@param x Name of column to check for shared substrings
#'(character(1), default: "value")
#'@param id Name of ID column uniquely identifying rows
#'(character(1), default: "ID")
#'@param new_col Name of column to be added to df
#'(character(1), default: "AB_group")
#'@importFrom dplyr arrange
#'@export
#'@returns Input data.frame df, with an new numeric column
#'(default name "AB_group"), where numbers correspond to groups.
#'@examples
#'df <- data.frame(Antigen=c("CD274 (PD-L1)", "PD-L1(B7-H1)", "B7-H1"), ID=1:3)
#'
#'# Format for using sharedSubstr
#'splitUnnest(df)
#'sharedSubstr(df, x="Antigen")
sharedSubstr <- function(df, x="value", id="ID", new_col="AB_group"){
    if (! all(c(x, id) %in% colnames(df))) {
        stop("x and id must be columns in df")
    }
    df <- dplyr::arrange(df, !!sym(id)) %>%
        dplyr::mutate(!!new_col := .sharedSubstr(df[[x]], df[[id]]))
    return(df)
}


# .sharedSubstr ----
#
#' Group based on shared substrings
#'
#' Takes a vector of words and a corresponding vector of group ids
#' and returns a numeric vector indicating groupings where members of a
#' group share a word with at least one other member of the group.
#'
#'@param x Character vector of words
#'@param id Vector of group ids
#'@returns A numeric vector the same length as the input vector x, where
#'each number is a group and groups are defined by sharing at least one
#'substring.
#'@examples
#'# "fox" occurs in groups 1 and 3, and "box" occurs in groups 3 and 4,
#'# so these groups are grouped together
#'AbNames:::.sharedSubstr(c("fox", "cat", "fox", "in", "box", "box"),
#'c(1, 2, 3, 3, 3, 4))
.sharedSubstr <- function(x, id){
    # NOTE: need to arrange query_df by ID so that 1 ID forms only one run

    # Make a co-occurrence matrix
    x_ln <- length(x)
    xx <- rep(x, x_ln)
    yy <- rep(x, each=x_ln)

    # Co-occurrence matrix, one row/column per ID value
    eq <- matrix(as.numeric(xx == yy), nrow = x_ln)
    eq <- rowsum(t(rowsum(eq, id, reorder=FALSE)), id, reorder=FALSE)

    # Find matches between groups
    gp_rle <- rle(id)
    gp <- seq_along(gp_rle$values)

    for (i in seq_along(gp_rle$values)){
        gp_eq <- gp[eq[i, ] > 0]
        gp[gp_eq] <- min(gp_eq)
    }

    # Convert back to (numeric) id vector
    result <- as.numeric(as.factor(rep(gp, gp_rle$lengths)))
    return(result)
}
