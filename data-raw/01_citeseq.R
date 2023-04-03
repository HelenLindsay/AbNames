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
citeseq_q <- citeseq %>% dplyr::select(ID, Antigen, Isotype_Control)
# Note - this does not search for control antibodies
query_df <- AbNames::makeQueryTable(citeseq_q,
                                    ab = "Antigen",
                                    control_col = "Isotype_Control")

# To allow re-running, remove ID columns before re-annotating -----

id_cols <- c("ALT_ID", "HGNC_ID", "HGNC_SYMBOL", "ENSEMBL_ID",
             "UNIPROT_ID", "ENTREZ_ID", "SOURCE")

citeseq <- citeseq %>%
    dplyr::select(-any_of(id_cols))

# Annotate using the gene aliases table ----

alias_results <- searchAliases(query_df)

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

# Check for missing annotation ----
citeseq %>%
    dplyr::filter(is.na(ALT_ID) & ! Isotype_Control)

# Patch the antigens that are still missing -----
cs_patch <- tibble::tribble(~Antigen, ~value,
                            "DopamineD4receptor", "dopamine receptor D4",
                            "DopamineReceptorD4", "dopamine receptor D4")
                            #"c-Fos", "FOS")
cs_patch <- cs_patch %>%
    dplyr::left_join(gene_aliases, by = "value") %>%
    dplyr::select(any_of(colnames(citeseq)))

citeseq <- citeseq %>%
    dplyr::rows_patch(cs_patch) %>%
    dplyr::ungroup()

# Fill Cat_Number, Oligo_ID, etc ----

# First fill other entries given (TotalSeq) catalogue number:
# If catalogue number is available, can fill clone, oligo
# For BD Bioscieces, I'm not sure if Antigen and Clone is enough to fill
# Cat_Number

exp_nrow <- nrow(citeseq)
n_na <- citeseq %>%
    dplyr::summarise(across(c("Cat_Number", "Clone", "Oligo_ID"),
                            ~sum(is.na(.x))))

citeseq_custom <- citeseq %>%
    dplyr::filter(grepl("custom", Cat_Number) | isTRUE(Custom_Antibody))


# TO DO: CHECK FOR ERRORS AND REMOVE "IGNORE"
citeseq <- citeseq %>%
    dplyr::filter(! grepl("custom", Cat_Number) &
                      (is.na(Custom_Antibody) | ! Custom_Antibody)) %>%
    fillByGroup(group = "Cat_Number",
                fill = c("Clone", "Oligo_ID", "TotalSeq_Cat"),
                multiple = "ignore") %>%
    # Allow leniency for Reactivity as it hasn't been sorted
    fillByGroup(group = "Cat_Number", fill = c("Reactivity"),
                multiple = "ignore") %>%
    fillByGroup(group = "RRID", fill = c("Clone", "Reactivity"),
                multiple = "ignore") %>%
    # It is possible to share oligo, antigen and clone but not cat number
    fillByGroup(group = c("Clone", "Oligo_ID", "TotalSeq_Cat"),
                fill = c("Cat_Number"), multiple = "ignore")

citeseq <- dplyr::bind_rows(citeseq, citeseq_custom) %>%
    dplyr::arrange(Study, Antigen)

exp_nrow == nrow(citeseq)

n_na <- citeseq %>%
    dplyr::summarise(across(c("Cat_Number", "Clone", "Oligo_ID"),
                            ~sum(is.na(.x))))


# Hao_2021 and Liu_2021 have entries for some but not all Cat_Numbers
# It appears info can be added for Liu, but Hao entries are custom
# (But do not ever have other group members)


# Standardising names for controls
#dplyr::group_by(Cat_Number) %>%
#dplyr::mutate(n = ifelse(n == 1 & grepl("kappa", Suggested_Antigen),
#                         100, n)) %>% # Prefer kappa if it's there


# Create citeseq data set ----
citeseq <- as.data.frame(citeseq)
usethis::use_data(citeseq, overwrite = TRUE, compress = "bzip2")


# HGNC:12102 = TRAV1-2 = TCR Va7.2?
# TRAV24, TRAJ18 = TCRVa24-Ja18
# HLA-DR = HLA-DRA?

# For antigens like TCR alpha/beta, KIR2DL1/S1/S3/S5, want to split like subunit
# CD66a_c_e
# CD18 associates with CD11a (LFA-1) CD11-b (Mac-1)
# CD3 (CD3E)__Nathan_2021 should match CD3D,E and G?

# filling in suggested antigen

#If entry above and below are the same, add the suggested antigen

#all_clones <- all_clones %>%
#    dplyr::ungroup() %>%
#    dplyr::mutate(before = lead(Suggested_Antigen),
#                  after = lag(Suggested_Antigen),
#                  Suggested_Antigen =
#                      ifelse(is.na(Suggested_Antigen) & before == after, before,
#                             Suggested_Antigen)) %>%
#    dplyr::select(-n, -before, -after, -AG_TEMP) %>%
#    dplyr::relocate(Antigen, Suggested_Antigen, Cat_Number, Clone)


# Cases:
# The antigen only appears in one study
# The antigen appears twice but has completely different names
# (choose the gene? Or keep the form that is "^CD[0-9]+$")
# CD314 - AN ERROR HAS PROPAGATED! FIXED FOR KOTLIAROV, WHICH STUDIES?
# If clone and Antigen name are identical but totalseq cat different, fill

# Notes ----
# Granja has an anti-mouse antibody but reactivity is Human?
#x %>% filter(Experiment == "Granja_2019", Antigen == "CD3")

