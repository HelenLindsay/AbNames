# filter_by_union ----
#'@title Select rows matching any column present in another data.frame
#'
#'@description Select values from a data.frame df matching any
#'column from another data.frame or a selection of row indices.  If a second
#'data frame and a selection of rows is provided, values from df matching
#'any value in df2[rows, ] are returned.  If only rows indices are provided,
#'rows matching any value in df[rows, ] are returned.  NAs are not matched.
#'
#'This is a wrapper around dplyr::filter with if_any for the case where the
#'value should be in a reference set, with the option to choose whether to use
#'another data frame or a subset of rows as the reference.
#'
#'@param df A data.frame from which to select matching rows
#'@param df2 Optional, a second data frame
#'@param rows Row indices for subsetting, either df2 if present or df
#'@param by = columns to select from df2
#'@author Helen Lindsay
#'@export
filter_by_union <- function(df, df2 = NULL, rows = NULL, by = NULL){
    if (! is.null(df2) & ! is.null(rows)){
        message("Row selection will be made from df2")
    }

    tmp <- .tempColName(df)
    qdf <- df
    if (! is.null(df2)) { qdf <- df2 }
    if (! is.null(rows)) qdf <- qdf[rows, ]

    # Check that values of by exist in qdf
    if (! is.null(by) & ! all(by %in% colnames(qdf))){
        stop("Not all columns in 'by' appear in df2:",
             toString(setdiff(by, colnames(qdf))))
    }

    # If by is not specified, use all shared columns
    if (is.null(by)) by <- intersect(colnames(qdf), colnames(df))

    # If there are names, assume they use dplyr join syntax, i.e. names
    # refer to columns in df, values to columns in df2
    if (! is.null(names(by))){
        no_name <- names(by) == ""
        names(by)[no_name] <- by[no_name]

        # Check that names of "by" exist in df
        if (! all(names(by) %in% colnames(df))){
            stop("Not all columns in 'by' appear in df:",
                 toString(setdiff(names(by), colnames(df))))
        }

        # Rename columns in qdf to match "by", set by to equal names of by
        qdf <- qdf[, by]
        colnames(qdf) <- names(by)
        by <- names(by)
    }

    keep_rows <- rowSums(vapply(by, function(x){
        df[[x]] %in% qdf[[x]] & ! is.na(df[[x]])
    }, logical(nrow(df)))) >= 1

    return(df[keep_rows,, drop = FALSE])
}


# group_by_any ----
#
# Group data.frame by a match in any of the grouping columns
#
# Not tested on more than 2 groups
# Pairwise comparison necessary?
#
# df - data.frame
# groups - character vector of grouping columns
# ignore e.g. Cat_Number == "custom_made".  Names should be column names
# in df, if they are not they are ignored. Regex?
#'@importFrom stats na.omit
group_by_any <- function(df, groups, new_col = "group", ignore = NULL,
                         verbose = FALSE){
    if (length(groups) < 2){
        if (isTRUE(verbose)){
            message("With only one group, group_by_any is ",
                    "equivalent to group_by\n")
        }
        result <- df %>%
            dplyr::group_by(!!rlang::sym(groups)) %>%
            dplyr::mutate(!!new_col := dplyr::cur_group_id())
        return(result)
    }

    idx <- purrr::map(groups,
                      ~dplyr::group_by(df, !!sym(.x)) %>%
                          dplyr::mutate(idx = dplyr::cur_group_id()) %>%
                          dplyr::pull(idx))

    idx <- do.call(dplyr::bind_cols, structure(idx, names = groups))
    idx[is.na(df[, groups])] <- NA

    # Remove groupings that should be ignored
    idx <- .ignore_groups(df, ignore, idx)

    new_idxs <- idx[[ groups[1] ]]
    curr_idxs <- idx[[ groups[1] ]]
    for (gp in groups[2:length(groups)]){
        merge_idxs <- idx[[gp]]
        for (val in sort(unique(na.omit(new_idxs)))){
            merge_vals <- na.omit(unique(merge_idxs[curr_idxs == val]))
            v_to_replace <- new_idxs[merge_idxs %in% merge_vals]
            if (any(v_to_replace)) {
                rp_id <- merge_idxs %in% merge_vals |
                    new_idxs %in% na.omit(v_to_replace)
                new_idxs[rp_id] <- min(v_to_replace, na.rm = TRUE)
            }
        }
    }

    df[[new_col]] <- new_idxs
    df %>% dplyr::group_by(!!sym(new_col))
}


# .ignore_groups ----
# Remove indices that should be ignored.  If names in "ignore" don't match
# columns of df, they are ignored
.ignore_groups <- function(df, ignore, idx){
    if (is.null(idx)) return(idx)

    for (nm in intersect(names(ignore), colnames(df))){
        vals <- ignore[[nm]]
        vals <- do.call(paste, list(vals, collapse = "|"))
        idx[grepl(vals, df[[nm]]), nm] <- NA
    }
    return(idx)
}

# left_join_any ----
# Join by matches in any set of columns
# Perform successive inner joins
# Allows cols in groups, e.g. "Cat_Number", c("Antigen", "Clone"), ...
# Cols must be a list.  Assume that rows can be uniquely identified
# shared: which function to use for shared columns?  rows_patch (fill NA) or
# rows_update (overwrite)
# ignore option for patching?  Set patch_cn to NULL?
# What to do if Antigen matches one row and Clone matches a different one?
# Probably should include both rows
# Do we expect columns to be NA if they do not match?  Not necessarily,
# Antigen may differ
# As is, only one column is given the chance to match
#'@author Helen Lindsay
#'@importFrom dplyr if_all
#'@importFrom dplyr everything
left_join_any <- function(x, y, cols, shared = c("patch", "update")){

    update_fun <- match.arg(shared)
    update_fun <- if (update_fun == "patch"){
        dplyr::rows_patch
    } else {
        dplyr::rows_update
    }

    # For each join, we only want to include the join column of interest in y
    cn_y <- colnames(y)
    cn_x <- colnames(x)
    join_cns <- unique(unlist(cols))

    if (! all(join_cns %in% cn_x) & all(join_cns %in% cn_y)){
        stop("Columns for joining must appear in both x and y")
    }

    # Columns to add are columns not already in in x and not used for joining
    # Shared columns are treated differently, by patching NA values
    patch_cn <- intersect(cn_x, setdiff(cn_y, join_cns))
    add_cn <- setdiff(cn_y, c(join_cns, patch_cn))
    aj_cn <- c()

    res <- vector("list", length = length(cols))

    # Inner join separately with each col set
    for (i in seq_along(cols)){
        col_set <- cols[[i]]

        y_subs <- y %>%
            dplyr::select(all_of(c(add_cn, col_set))) %>%
            unique()

        new_res <- x %>%
            dplyr::inner_join(y_subs, by = col_set,
                              na_matches = "never", multiple = "all",
                              relationship = "many-to-many")

       # If columns are shared, either update or patch values in x from y
        if (length(patch_cn) > 0){
            y_patch <- y %>%
                dplyr::select(all_of(c(patch_cn, col_set))) %>%
                dplyr::filter(if_all(everything(), complete.cases)) %>%
                unique()
            new_res <- update_fun(new_res, y_patch,
                                  unmatched = "ignore", by = col_set)
        }

        res[[i]] <- new_res
    }

    res <- do.call(dplyr::bind_rows, res)
    x_aj <- dplyr::anti_join(x, res, by = setdiff(cn_x, patch_cn))

    # Add the unmatched rows back in
    res <- dplyr::bind_rows(x_aj, res) %>%
        unique()

    return(res)

}


# rows_patch_any ----
#
# simple wrapper for coalesce, rearranging rows
# match returns the first match, there may be many
# only works for a single match column
# apply where, to match columns, coalesce, discard if multiple values?
rows_patch_any <- function(x, y, match_col){
    row_order <- match(x[[match_col]], y[[match_col]])
    y <- y[row_order,]

}



