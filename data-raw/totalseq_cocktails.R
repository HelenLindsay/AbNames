# Setup ----

library(readxl)
library(dplyr)
library(readr)
library(AbNames)

# Download TotalSeq cocktail information ----

bl <- "https://www.biolegend.com/Files/Images/BioLegend/totalseq/"

tsa <- paste0(bl, "TotalSeq_A_Human_Universal_Cocktail_v1_163_",
             "Antibodies_399907_Barcodes.xlsx")

tsb <- paste0(bl, "TotalSeq_B_Universal_Cocktail_v1_140_",
              "Antibodies_399904_Barcodes.xlsx")

tsc <- paste0(bl, "TotalSeq_C_Human_Universal_Cocktail_v1_137_",
              "Antibodies_399905_Barcodes.xlsx")

tsd <- paste0(bl, "TotalSeq_D_Human_Heme_Oncology_Cocktail_V01_Barcode.xlsx")

totalseq_fnames <- c(tsa, tsb, tsc, tsd)

totalseq <- lapply(totalseq_fnames, function(fn){
    f <- tempfile()
    download.file(fn, destfile = f)
    readxl::read_xlsx(f)
})


nrows <- sapply(totalseq, nrow)

# Get the TotalSeq category from the filename and add to tables ----
ts_cat <- gsub(".*TotalSeq_([ABCD])_.*", "\\1", basename(totalseq_fnames))
ts_cat <- rep(ts_cat, nrows)
ts_cat <- split(ts_cat, ts_cat)

totalseq <- lapply(seq_along(totalseq), function(i) {
    totalseq[[i]] %>% dplyr::mutate(TotalSeq_Cat = ts_cat[[i]])
})

totalseq <- Reduce(dplyr::full_join, totalseq)

# Check that no rows have been lost or gained ----
nrow(totalseq) == sum(nrows)

# Coalesce columns with the same meaning ----
totalseq <- totalseq %>%
    dplyr::rename(Ensembl_ID = `Ensembl ID`,
                  Gene_Symbol = `Gene Name`,
                  Barcode_Sequence = `Barcode sequence`,
                  Oligo_ID = DNA_ID,
                  Antigen = Description) %>%
    dplyr::mutate(Barcode_Sequence =
                      dplyr::coalesce(Barcode_Sequence, Barcode, Sequence),
                  Oligo_ID = dplyr::coalesce(Oligo_ID, `Format / Barcode`),
                  Clone = dplyr::coalesce(Clone, clone),
                  Ensembl_ID = dplyr::coalesce(Ensembl_ID, `Ensemble ID`),
                  Gene_Symbol = dplyr::coalesce(Gene_Symbol, `Gene name`),
                  Antigen = dplyr::coalesce(Antigen, Specificity)) %>%
    dplyr::select(-Barcode, -`Format / Barcode`, -clone,
                  -`Ensemble ID`, -`Gene name`, -Sequence, -Specificity) %>%
    dplyr::mutate(Reactivity = ifelse(grepl("anti-human\\s", Antigen),
                                      "human", NA),
                  Antigen = gsub("anti-human\\s", "", Antigen),
                  Oligo_ID = substr(Oligo_ID, 2, nchar(Oligo_ID)),
                  Antigen = AbNames::replaceGreekSyms(Antigen, "sym2letter"))%>%
    dplyr::relocate(Antigen, Clone, Ensembl_ID, Gene_Symbol, Oligo_ID,
                    TotalSeq_Cat, Barcode_Sequence, Reactivity)

# Some TotalSeq B Ensembl_IDs are duplicated barcode sequences, set to NA ----
totalseq <- totalseq %>%
    dplyr::mutate(Ensembl_ID = ifelse(grepl("^ENSG", Ensembl_ID),
                                      Ensembl_ID, NA))

# Check that there is only one Ensembl_ID per Antigen-Clone combo ----

totalseq <- totalseq %>%
    AbNames::nPerGroup(group = c("Antigen", "Clone"), "Ensembl_ID")

if (max(totalseq$n_per_group) > 1){
    warning("More than one value of Ensembl_ID per Antigen/Clone combo")
}

totalseq <- totalseq %>% dplyr::select(-n_per_group)

# Fill in missing genes if Antigen, Clone, and Oligo_ID match ----
# Note that by adding "Reactivity", controls will be removed as reactivity was
# only defined for human
totalseq <- totalseq %>%
    AbNames::fillByGroup(group = c("Antigen", "Clone", "Reactivity"),
                         fill = c("Ensembl_ID", "Gene_Symbol"))

# Fill in the missing gene symbols if Ensembl ID is given
totalseq <- totalseq %>%
    AbNames::fillByGroup(group = "Ensembl_ID", fill = "Gene_Symbol")


# Create totalseq_cocktails data set ----
totalseq_cocktails <- as.data.frame(totalseq)
usethis::use_data(totalseq_cocktails, overwrite = TRUE, compress = "bzip2")
