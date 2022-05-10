# searchHGNC ----
#'@title searchHGNC
#'@description Search HGNC gene symbols, aliases and names for an exact match
#'to an value, assumed to be an antigen name or part of an antigen name.
#'@param df A data.frame (or similar object, e.g. a tibble)
#'@importFrom utils data
searchHGNC <- function(df){
    utils::data("hgnc_long", envir = environment())

}
