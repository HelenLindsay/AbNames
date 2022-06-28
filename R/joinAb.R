# join_ab -----
#
# Join two data.frames by corresponding but non-identical antibody columns
#
# Uses a combination of string distance and common "words".  Intended for
# matching vendor antibody names with names provided in a data set.

join_ab <- function(x, y, x_col = "Antigen", y_col = "Data_Name"){
    rlang::check_installed("fuzzyjoin")

    xtmp <- .tempColName(x)
    ytmp <- .tempColName(y)
    x <- dplyr::mutate(x, !!xtmp := replacePunct(.data[[x_col]]))
    y <- dplyr::mutate(y, !!ytmp := replacePunct(.data[[y_col]]))
    by = setNames(ytmp, xtmp)
    x_join <- strdist_join(x, y, by = by)

}


replacePunct <- function(x){
    return(gsub("[[:punct:]]", "-", x))
}


# Uses fuzzyjoin::stringdist_left_join, which returns all matches within
# a given max_dist
strdist_join <- function(x, y, by, max_dist = 2){
    sdist_col <- .tempColName(x, 1, "strdist")
    x <- x %>%
        fuzzyjoin::stringdist_left_join(y, max_dist = max_dist, by = by,
                                         distance_col = sdist_col) %>%
        dplyr::rename()


         dplyr::rename(Antigen = `Antigen.x`) #%>%
    #     dplyr::select(-`Antigen.y`) %>%
    #     dplyr::relocate(Data_Name) %>%
    #     dplyr::group_by(Antigen) %>%
    #     dplyr::filter(strdist == min(strdist)) %>%
    #     dplyr::group_by(Data_Name) %>%
    #     dplyr::filter(strdist == min(strdist)) %>%
    #     dplyr::filter(! is.na(Data_Name)) %>%
    #     dplyr::select(-strdist)
    #
    # cn <- colnames(triana)
    # triana <- triana %>%
    #     dplyr::rows_patch(triana_um, by = setdiff(cn, "Data_Name"))
         # dplyr::rows_patch(triana_um, by = setdiff(cn, "Data_Name"))
         # dplyr::rows_patch(triana_um, by = setdiff(cn, "Data_Name"))
         # dplyr::rows_patch(triana_um, by = setdiff(cn, "Data_Name"))
}
