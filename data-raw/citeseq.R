library("tidyverse")
library("AbNames")
library("janitor")
library("stringi")

# To do: make sure human IDs are not assigned to non-human reactive Abs
# To do: group and check differences in ALT_ID.
# Chung c-Fos has wrong catalogue number?

citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
citeseq <- read.csv(citeseq_fname) %>% unique()

# Remove non-ascii characters -----

citeseq <- citeseq %>%
    dplyr::mutate(across(where(is.character),
                         ~stringi::stri_trans_general(.x,
                                id="Any-Latin;Greek-Latin;Latin-ASCII")))


citeseq <- AbNames::addID(citeseq)
citeseq_q <- citeseq %>% dplyr::select(ID, Antigen)
query_df <- AbNames::makeQueryTable(citeseq_q, ab = "Antigen")

alias_results <- searchAliases(query_df)

id_cols <- c("ALT_ID", "HGNC_ID", "HGNC_SYMBOL", "ENSEMBL_ID",
             "UNIPROT_ID", "ENTREZ_ID", "SOURCE")

# Remove matches to several genes, select just columns of interest

alias_results <- alias_results %>%
    # Select ID and HGNC columns
    dplyr::select(matches("ID|HGNC"), name, SOURCE) %>%
    unique() %>% # Collapse results with same ID from different queries
    group_by(ID) %>%
    # Collapse multi-subunit entries, convert "NA" to NA
    dplyr::summarise(dplyr::across(all_of(id_cols), ~toString(unique(.x)))) %>%
    dplyr::mutate(dplyr::across(all_of(id_cols), ~na_if(.x, "NA")))

nrow(alias_results)

citeseq <- citeseq %>%
    dplyr::left_join(alias_results, by = "ID") %>%
    dplyr::relocate(ID, Antigen, Cat_Number, HGNC_ID) %>%
    unique()

citeseq <- searchTotalseq(citeseq)

# Patch the antigens that are still missing -----
cs_patch <- tibble::tribble(~Antigen, ~value,
                            "DopamineD4receptor", "dopamine receptor D4",
                            "DopamineReceptorD4", "dopamine receptor D4",
                            "c-Fos", "FOS")
cs_patch <- cs_patch %>%
    dplyr::left_join(gene_aliases, by = "value") %>%
    dplyr::select(any_of(colnames(citeseq)))

citeseq <- citeseq %>%
    dplyr::rows_patch(cs_patch)

# Create citeseq data set ----
citeseq <- as.data.frame(citeseq)
usethis::use_data(citeseq, overwrite = TRUE, compress = "bzip2")


# TCR alpha/beta - not matched correctly?
# HGNC:12102 = TRAV1-2 = TCR Va7.2?
# TRAV24, TRAJ18 = TCRVa24-Ja18
# HLA2 != HLA-A
# HLA-DR = HLA-DRA?
# HLA.A.B.C != HGNC:914

# For antigens like TCR alpha/beta, KIR2DL1/S1/S3/S5, want to split like subunit
# CD66a_c_e
# CD18 associates with CD11a (LFA-1) CD11-b (Mac-1)
# CD3 (CD3E)__Nathan_2021 should match CD3D,E and G?


