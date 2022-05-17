## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
library(AbNames)
library(dplyr)

data("hgnc", package = "AbNames") # hgnc_long is loaded similarly

# Show the first entries of hgnc, where each row is the start of one column
dplyr::glimpse(hgnc) 

# (Note: it isn't necessary to use dplyr:: to call "glimpse" as dplyr is loaded
# with the library call above. This syntax is used to make it clear which
# packages functions belong to)

## -----------------------------------------------------------------------------
data("totalseq_cocktails", package = "AbNames")
dplyr::glimpse(totalseq_cocktails)

## ---- message = FALSE---------------------------------------------------------
# As the citeseq data set contains raw data, it is loaded differently
# than the other data sets 

citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
citeseq <- read.csv(citeseq_fname)
dplyr::glimpse(citeseq)

## -----------------------------------------------------------------------------
# The regular expression inside "grepl" searches for PD.L1, PD-L1, PDL1 or CD274

cd274 <- citeseq %>%
    dplyr::filter(grepl("PD[\\.-]?L1|CD274", Antigen)) %>%
    dplyr::pull(Antigen) %>%
    unique()

cd274

## -----------------------------------------------------------------------------
# Add an ID column to the citeseq data set to allow the results table to be merged
citeseq <- AbNames::addID(citeseq)

# Select just the columns that are needed for querying:
citeseq_q <- citeseq %>% dplyr::select(ID, Antigen)

# Apply the default transformation sequence and make a query table in long format
query_df <- AbNames::makeQueryTable(citeseq_q, ab = "Antigen")

# Print one example antigen from the query table:
query_df %>%
    dplyr::filter(ID == "TCR alpha/beta-Hao_2021")


## -----------------------------------------------------------------------------
# Get the default sequence of formatting functions 
default_funs <- AbNames::defaultQuery()

# Print the first two formatting functions as an example
default_funs[1:2]

## -----------------------------------------------------------------------------
tcr <- data.frame(Antigen =
                      c("TCR alpha/beta", "TCRab", "TCR gamma/delta", "TCRgd",
                        "TCR g/d", "TCR Vgamma9", "TCR Vg9", "TCR Vd2",
                        "TCR Vdelta2", "TCR Vα24-Jα18", "TCRVa24.Ja18", 
                        "TCR Valpha24-Jalpha18", "TCR Vα7.2", "TCR Va7.2",
                        "TCRa7.2", "TCRVa7.2", "TCR Vbeta13.1", "TCR γ/δ", 
                        "TCR Vβ13.1", "TCR Vγ9", "TCR Vδ2", "TCR α/β",
                        "TCRb", "TCRg"))


# First, we convert the Greek symbols to letters.  
# Note that as we are using replaceGreekSyms in a dplyr pipeline, we don't 
# put quotes around the column name, i.e. Antigen not "Antigen".

tcr <- tcr %>%
    dplyr::mutate(query = 
                    AbNames::replaceGreekSyms(Antigen, replace = "sym2letter"))
 
# Print out a few rows to see the result
tcr %>% 
    dplyr::filter(Antigen %in% c("TCR Vβ13.1", "TCR Vδ2", "TCR α/β"))


## -----------------------------------------------------------------------------
tcr_f <- AbNames::formatTCR(tcr, tcr = "query")

# Print out the first few rows
tcr_f %>%
    head() 


## ---- echo = FALSE------------------------------------------------------------
set.seed(549867)

## -----------------------------------------------------------------------------
hgnc_results <- searchHGNC(query_df)

# Print 10 random results:
hgnc_results %>%
    dplyr::select(ID, name, value, symbol_type) %>%
    dplyr::ungroup() %>%
    dplyr::sample_n(10)

## -----------------------------------------------------------------------------
hgnc_results <- hgnc_results %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(n_ids = dplyr::n_distinct(HGNC_ID)) # Count distinct IDs
    
# Select antigens where there are matches to multiple genes but not
# because the antibody is against a multi-gene protein 
multi_gene <- hgnc_results %>%
    dplyr::filter(n_ids > 1, ! all(name %in% c("TCR_long", "subunit"))) %>%
    dplyr::select(ID, name, value, HGNC_ID)

# Look at the first group in multi_gene.
# Set interactive = FALSE for interactive exploration
showGroups(multi_gene, 1, interactive = FALSE)

# This is an example where the antibody is against a heterodimeric protein.
# We can confirm this by looking up the vendor catalogue number:
citeseq %>%
    dplyr::filter(ID == "CD11a.CD18-Wu_2021_b") %>%
    dplyr::select(Antigen, Cat_Number, Vendor)


## -----------------------------------------------------------------------------
hgnc_results %>%
    dplyr::filter(ID == "CD270 (HVEM, TR2)-Hao_2021") %>%
    data.frame()

## -----------------------------------------------------------------------------
nrow(hgnc_results)

# Remove matches to several genes, select just columns of interest
hgnc_results <- hgnc_results %>%
    dplyr::filter(n_ids == 1 | ! all(name %in% c("TCR_long", "subunit"))) %>%
    dplyr::select(ID, HGNC_ID, ENSEMBL_ID, UNIPROT_IDS)

nrow(hgnc_results)

citeseq <- citeseq %>%
    dplyr::full_join(hgnc_results, by = "ID") %>%
    dplyr::relocate(ID, Antigen, Cat_Number, HGNC_ID) %>%
    unique()

head(citeseq)

## -----------------------------------------------------------------------------
# Select some data to demonstrate filling:
# Get entries sharing the same catalogue number,
# where not every entry has a match
fill_demo <- citeseq %>%
    dplyr::group_by(Cat_Number) %>%
    dplyr::arrange(Cat_Number) %>%
    dplyr::filter(!is.na(Cat_Number), 
                  any(is.na(HGNC_ID)),
                  ! all(is.na(HGNC_ID)))

AbNames::showGroups(fill_demo, interactive = FALSE)

## -----------------------------------------------------------------------------
fill_demo <- AbNames::fillByGroup(fill_demo, "Cat_Number",
                                  fill = c("HGNC_ID", "TotalSeq_Cat", "Vendor",
                                           "ENSEMBL_ID", "UNIPROT_IDS"))
# Print out the first group again
AbNames::showGroups(fill_demo, interactive = FALSE)

## -----------------------------------------------------------------------------
sessionInfo()

