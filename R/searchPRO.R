# searchPRO ----
#
# Search the protein ontology for a match to an antibody
#
# As this function slow, we recommend only using this function for
# antibodies that were not matched by other means.
searchPRO <- function(df, ab, group = "ID", interactive = FALSE,
                      species = "human"){
    tmp <- .tempColName(df)

    utils::data("pro_long", envir = environment())

    exact_match <- df[[ab]] %in% pro_long$value
    sp_match <- sprintf("%s (%s)", df[[ab]], species) %in% pro_long$value

    df <- dplyr::mutate(df, !!tmp := exact_match | sp_match) %>%
        dplyr::group_by(!!rlang::sym(group))

    has_match <- dplyr::filter(df, any(!!rlang::sym(tmp)))
    df <- dplyr::filter(df, ! any(!! rlang::sym(tmp)))

    df_ab <- unique(df[[ab]])

}
