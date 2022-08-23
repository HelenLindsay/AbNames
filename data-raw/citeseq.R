library("tidyverse")
library("AbNames")

citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
citeseq <- read.csv(citeseq_fname) %>% unique()
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

# Create citeseq data set ----
citeseq <- as.data.frame(citeseq)
usethis::use_data(citeseq, overwrite = TRUE, compress = "bzip2")

# TCR alpha/beta - not matched correctly
# HGNC:12102 = TRAV1-2 = TCR Va7.2?
# TRAV24, TRAJ18 = TCRVa24-Ja18
# HLA2 != HLA-A
# HLA-DR = HLA-DRA?
# HLA.A.B.C != HGNC:914

split_merge_str <- sprintf('grepl("TCR", %s)', ab)

qdf <- citeseq %>% select(ID, Antigen) %>% unique()

qry = list(purrr::partial(gsubAb, ab = !!ab), # Remove A/antis
           purrr::partial(gsubAb, ab = !!ab, pattern = "\\s[Rr]ecombinant"),
           purrr::partial(splitUnnest, ab = !!ab), # Brackets
           purrr::partial(splitUnnest, ab = !!ab, split = ", "),  # Commas
           # / _ or . if at least 3 on each side and not TCR
           purrr::partial(splitUnnest, exclude = "TCR",
                          split = "(?<=[A-z0-9-]{3})[\\/_\\.](?=[A-z0-9-]{3,})"),
           purrr::partial(.reformatAb, ab = !!ab),
           purrr::partial(.reformatAb, ab = !!ab),
           purrr::partial(separateSubunits, ab = !!ab, new_col = "subunit"),
           purrr::partial(splitMerge, ex = !!split_merge_str,
                          f = formatTCR, tcr = "greek_letter")
)

qdf <- magrittr::freduce(qdf, qry)

# For antigens like TCR alpha/beta, KIR2DL1/S1/S3/S5, want to split like subunit
# CD66a_c_e
# CD18 associates with CD11a (LFA-1) CD11-b (Mac-1)
# CD3.2__Mimitou_2021_cite_and_dogma_seq in query table, only matches through
# totalseq = CD3E
#  CD3 (CD3E)__Nathan_2021 should match CD3D,E and G?
# Pombo DopamineD4receptor = "dopamine receptor D4" = "HGNC:3025


