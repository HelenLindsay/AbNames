# .warnIfColExists ----
#' Warn if a column already exists in a data.frame
#'
#'@param df A data.frame or tibble
#'@param new_col character(1) Name of a column to add to df
.warnIfColExists <- function(df, new_col){
    # Check if new_col already exists in df, stop if so

    if (new_col %in% colnames(df)){
        msg1 <- sprintf("Column %s already exists in data.frame.\n", new_col)
        msg2 <- "Please supply a different value for new_col"
        stop(msg1, msg2)
    }
}

# .tempColName ----
#' Create a name for a temporary column in a data.frame
#'
#' Returns a name that does not already exist in columns of a data.frame
#'@param df A data.frame or tibble
#'@param n (default 1) How many temporary column names are needed?
.tempColName <- function(df, n = 1){
    cn <- colnames(df)
    temp_col <- "TEMP"
    i <- 1
    res <- character()

    while(length(res) < n){
        if (! temp_col %in% cn) res <- c(res, temp_col)
        temp_col <- sprintf("TEMP%s", i)
        i <- i + 1
    }

    return(res)
}



