## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message = FALSE, warning=FALSE------------------------------------------

library("scRNAseq")
library("AbNames")

kotliarov <- KotliarovPBMCData(mode = "adt")
head(rownames(kotliarov))

## -----------------------------------------------------------------------------
# Create a data.frame for matching to the CITEseq data set
# Column "Antigen" will be used for matching
df <- data.frame(Original = rownames(kotliarov),
                 # Remove _PROT (optionally with a space before)
                 Antigen = gsub(" ?_PROT", "", rownames(kotliarov)))
head(df)


## -----------------------------------------------------------------------------
df <- matchToCiteseq(df)
print(head(df))

## -----------------------------------------------------------------------------
df[! df$Antigen == df$Antigen_std, ]

## -----------------------------------------------------------------------------
head(df[order(df$n_matched),], 10)

## -----------------------------------------------------------------------------

# To rename the singleCellExperiment, we need a named vector of the new names,
# with names being the original names of the singleCellExperiment.

# "structure" is a method that allows us to create a named vector in one step
new_nms <- structure(df[[2]], names = rownames(kotliarov))
print(head(new_nms))

kotliarov <- renameADT(kotliarov, new_nms)
head(rownames(kotliarov))


## ---- echo = FALSE------------------------------------------------------------
# Load CiteFuse example data

library(CiteFuse)
data("CITEseq_example", package = "CiteFuse")
sce_citeseq <- preprocessing(CITEseq_example)
rownames(altExp(sce_citeseq, "ADT"))

# DO WE NEED TO DO ANY PREPROCESSING OF THE NAMES?  JUST PASS IN CITESEQ?
df <- data.frame(Antigen = rownames(altExp(sce_citeseq, "ADT")))
df <- matchToCiteseq(df)

# Print entries that differ from the original
df[! df$Antigen == df$Antigen_std, ]

# Print the entries with the fewest matches
head(df[order(df$n_matched),], 10)


## ---- message = FALSE---------------------------------------------------------
# We load dplyr to simplify data.frame manipulations.
library(dplyr)

# Load citeseq data set
data(citeseq)

tcrb <- citeseq %>%
    # We will first select entries that equal "TCRb"
    dplyr::filter(Antigen == "TCRb") %>%
    # Then subset to just the columns we want to use for matching
    dplyr::select(Antigen, Clone, Cat_Number, ALT_ID)

# Then we will find entries in citeseq that match any of these columns.

# filter_by_union is a function to return rows where any the value of any
# column occurs in a reference data.frame
citeseq %>%
    filter_by_union(tcrb) %>%
    # Select just columns of interest for viewing
    dplyr::select(Study, Antigen, Clone, Cat_Number, ALT_ID)


## -----------------------------------------------------------------------------
sessionInfo()

## ---- echo = FALSE------------------------------------------------------------
# alternative data sets
#library(SingleCellMultiModal)
#cord_blood <- CITEseq(DataType="cord_blood", modes="scADT_Counts",
#                      dry.run=FALSE, DataClass = "SingleCellExperiment")
#cord_blood



