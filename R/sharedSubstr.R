# sharedSubstr ----
#
#' Group based on shared substrings
#'
#' Takes a vector of words and a corresponding vector of group ids
#' and returns a numeric vector indicating groupings where members of a
#' group share a word with at least one other member of the group.
#'
#' The intention is that if an antibody has alternative names and is sometimes
#' called by both, e.g. CD274 (B7-H1), we would like to group all entries
#' matching CD274 with all entries matching B7-H1.
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
sharedSubstr <- function(df, x = "value", id = "ID", new_col = "AB_group"){
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
    yy <- rep(x, each = x_ln)

    # Co-occurrence matrix, one row/column per ID value
    eq <- matrix(as.numeric(xx == yy), nrow = x_ln)
    eq <- rowsum(t(rowsum(eq, id, reorder = FALSE)), id, reorder = FALSE)

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


# #'@importFrom dplyr cur_group_id
#sharedSubstrDf <- function(df, id = "ID", x = "value"){
#    occ_df <- dplyr::group_by(df, !!sym(x)) %>%
#        dplyr::mutate(n = cur_group_id())
#
#    n_words <- max(occ_df$n)
#
#    occ_df <- dplyr::group_by(occ_df, !!sym(id)) %>%
#        dplyr::summarise(n = list(unique(n)))
#
#    # Matrix rows are IDs, columns are words
#    m <- matrix(0, nrow = nrow(occ_df), ncol = n_words)
#    word_occ <- pull(occ_df, n)
#
#    for (i in seq_len(nrow(occ_df))){
#        m[i,  word_occ[[i]]] <- 1
#    }
#
#    # Assign group IDs to IDs based on shared words
#    ids <- seq_len(nrow(occ_df))
#    for (i in seq_len()
#
#
#    (i in seq_len(n_words)){
#        m[i, ]
#    }

# Why no match for CD158f (KIR2DL5)__Su_2020 CD158f?
    # because KIR2DL5 matches A and B
# CD16?? = FCGR3A (but also FCGR3B??), cat number 302061
    # correction comes from totalseq
# Check for PD-1 / HGNC = PDCD1 - some filled by symbol?
# CHECK CD3, some are CD3D, some CD3E - 300475 filled to CD3D??
# - not filled for Stuart, no TotalSeq_Cat or vendor?

# Careful with CD3.2 - Mimitou, not equal to CD32!
# HLA-DR / HLA-DRA HGNC:4947 good example!



