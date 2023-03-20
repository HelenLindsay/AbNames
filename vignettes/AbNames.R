## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
library(AbNames)
library(dplyr)

data("gene_aliases", package = "AbNames") 

# Show the first entries of gene_aliases,
# where each row is the start of one column
dplyr::glimpse(gene_aliases) 

# (Note: it isn't necessary to use dplyr:: to call "glimpse" as dplyr is loaded
# with the library call above. This syntax is used to make it clear which
# packages functions belong to)

## -----------------------------------------------------------------------------
data("totalseq", package = "AbNames")
dplyr::glimpse(totalseq)

## ---- message = FALSE---------------------------------------------------------
# As the citeseq data set contains raw data, it is loaded differently
# than the other data sets 

citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
citeseq <- read.csv(citeseq_fname) %>% unique()
dplyr::glimpse(citeseq)

# TO DO: AS THIS IS ONLY A SUBSET OF THE CITESEQ INFO, THERE ARE APPARENT
# DUPLICATES, DECIDE WHAT TO DO WITH THESE

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

# Remove control columns, as we are only searching human proteins
controls <- dplyr::filter(citeseq, Control)
citeseq <- dplyr::filter(citeseq, ! Control)

# Select just the columns that are needed for querying:
citeseq_q <- citeseq %>% dplyr::select(ID, Antigen)

# Apply the default transformation sequence and make a query table in long format
query_df <- AbNames::makeQueryTable(citeseq_q, ab = "Antigen")

# Print one example antigen from the query table:
query_df %>%
    dplyr::filter(ID == "TCR alpha/beta__Hao_2021")


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

