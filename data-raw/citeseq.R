library("tidyverse")
library("AbNames")

citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
citeseq <- read.csv(citeseq_fname) %>% unique()
citeseq <- AbNames::addID(citeseq)
citeseq_q <- citeseq %>% dplyr::select(ID, Antigen)
query_df <- AbNames::makeQueryTable(citeseq_q, ab = "Antigen")
alias_results <- searchAliases(query_df)

id_cols <- c("HGNC_ID", "HGNC_SYMBOL", "ENSEMBL_ID", "UNIPROT_ID")

# Remove matches to several genes, select just columns of interest

alias_results %>%
    dplyr::select(matches("ID|HGNC"), name) %>% # Select ID and HGNC columns
    unique() %>% # Collapse results with same ID from different queries

    # Collapse multi-subunit entries, convert "NA" to NA
    dplyr::summarise(dplyr::across(all_of(id_cols), ~toString(unique(.x)))) %>%
    dplyr::mutate(dplyr::across(all_of(id_cols), ~na_if(.x, "NA")))

nrow(alias_results)

citeseq <- citeseq %>%
    dplyr::left_join(alias_results, by = "ID") %>%
    dplyr::relocate(ID, Antigen, Cat_Number, HGNC_ID) %>%
    unique()

citeseq <- searchTotalseq(citeseq)

# Create citeseq data set ----
citeseq <- as.data.frame(citeseq)
usethis::use_data(citeseq, overwrite = TRUE, compress = "bzip2")
