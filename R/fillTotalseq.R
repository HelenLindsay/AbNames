# fillTotalseq ----
#
#' Fill empty gene identifiers using TotalSeq information
#'
#'@param df A data.frame containing column(s) to fill with gene identifiers.
#'@param ab character(1) Name of the column in df containing antibody names.
#'@importFrom dplyr all_of
fillTotalseq <- function(df, ab = "Antigen"){
    utils::data("totalseq_cocktails", envir = environment())
    common_cols <- intersect(colnames(df), colnames(totalseq_cocktails))
    ts <- totalseq_cocktails %>%
        dplyr::select(all_of(common_cols))


}
