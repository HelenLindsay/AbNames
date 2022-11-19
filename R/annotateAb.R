# Annotate antigens with identifiers
#
# This is a wrapper for the default matching pipeline using functions
# makeQueryTable and searchAliases.
#
annotateAb <- function(x, control_col = NULL){
    x <- AbNames::addID(x)

    # Remove control columns, as we are only searching human proteins
    controls <- dplyr::filter(x, Control)
    x <- dplyr::filter(x, ! Control)

    # Select just the columns that are needed for querying:
    citeseq_q <- x %>% dplyr::select(ID, Antigen)

    # Apply the default transformation sequence and make a query table in long format
    query_df <- AbNames::makeQueryTable(x, ab = "Antigen")



    alias_results <- alias_results %>%
        dplyr::group_by(ID) %>%
        dplyr::mutate(n_ids = dplyr::n_distinct(HGNC_ID)) # Count distinct IDs

    # Select antigens where there are matches to multiple genes but not
    # because the antibody is against a multi-gene protein
    multi_gene <- alias_results %>%
        dplyr::filter(n_ids > 1, ! all(name %in% c("TCR_long", "subunit"))) %>%
        dplyr::select(ID, name, value, HGNC_ID)

    alias_results %>%
        dplyr::select(matches("ID|HGNC"), name) %>% # Select ID and HGNC columns
        unique() %>% # Collapse results with same ID from different queries

        # Collapse multi-subunit entries, convert "NA" to NA
        dplyr::summarise(dplyr::across(all_of(id_cols), ~toString(unique(.x)))) %>%
        dplyr::mutate(dplyr::across(all_of(id_cols), ~na_if(.x, "NA")))

    nrow(alias_results)

    x <- x %>%
        dplyr::left_join(alias_results, by = "ID") %>%
        unique()

    ts <- totalseq %>%
        dplyr::select(any_of(colnames(citeseq)))

    missing <- count_missing(citeseq)

    # Fill in IDs where Antigen, Oligo, Clone and TotalSeq category match
    x <- x %>%
        dplyr::rows_patch(ts,
                          by = c("Antigen", "Oligo_ID", "Clone", "TotalSeq_Cat"),
                          unmatched = "ignore")


    x <- x %>%
        dplyr::mutate(across(any_of("HGNC_ID", "ENSEMBL_ID",)))
}
