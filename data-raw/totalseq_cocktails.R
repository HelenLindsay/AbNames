# Notes ----

# TotalSeq_Barcodes is more complete than TotalSeq_Cocktails.
# TotalSeq_Cocktails just used to fill in missing isotype control names.

# Setup ----

library(readxl)
library(dplyr)
library(readr)
library(AbNames)
library(stringr)

# Download TotalSeq Barcode lookup information ----

# (Manually downloaded from
# https://www.biolegend.com/en-us/totalseq/barcode-lookup
# and converted to csv due to import error)

# Note: "NA" antigens are isotype controls and streptavidin conjugates

dest <- "~/Analyses/CITEseq_curation/data/totalseq_barcodes"
dest_fnames <- c("totalseq_a_antibodies.csv", "totalseq_b_antibodies.csv",
                 "totalseq_c_antibodies.csv", "totalseq_d_antibodies.csv")
dest_fnames <- file.path(dest, dest_fnames)

totalseq_cat <- c("A", "B", "C", "D")

ts_barcodes <- lapply(seq_along(dest_fnames), function(i){
    readr::read_delim(dest_fnames[i]) %>%
        dplyr::mutate(TotalSeq_Cat = totalseq_cat[i])
})

ts_barcodes <- do.call(bind_rows, ts_barcodes) %>%
    dplyr::rename(Barcode_Sequence = Sequence,
                  ENSEMBL_ID = `Ensembl Gene Id`,
                  Date_Released = `Date Released`,
                  Cat_Number = `Catalog`,
                  Oligo_ID = Barcode,
                  Antigen = Description) %>%
    dplyr::relocate(TotalSeq_Cat, .after = Clone) %>%
    tidyr::separate(Reactivity, into = c("Reactivity", "Cross_Reactivity"),
                    sep =
            ", (<[Bb]>|<strong>)?Cross-Reactivity:(</[Bb]>|</strong>)? ",
            fill = "right") %>%
    dplyr::mutate(across(where(is.character), ~ stringr::str_squish(.x)),
                  ENSEMBL_ID = gsub(";|:|nbsp", "", ENSEMBL_ID)) %>%
    dplyr::filter(! (is.na(Antigen) & is.na(Clone)))

# Check if one oligo is assigned to exactly one
# Antigen / Clone / TotalSeq_Cat comb

# (The only entry is same antigen written differently TCR Vα2 / TCR V&alpha;2)
x <- ts_barcodes %>%
    nPerGroup(group = "Oligo_ID",
              col = c("Antigen", "Clone", "TotalSeq_Cat")) %>%
    filter(if_any(c(nAntigen, nClone), ~.x > 1))

# From this information, it appears that Antigen / Oligo / TotalSeq_Cat is
# enough to fill in Cat_Number and Clone

# Is Clone enough to fill in Antigen?


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

totalseq <- Reduce(dplyr::full_join, totalseq) %>%
    dplyr::mutate(across(.cols = everything(), stringr::str_squish))

# Check that no rows have been lost or gained ----
nrow(totalseq) == sum(nrows)

# Coalesce columns with the same meaning ----
totalseq <- totalseq %>%
    dplyr::rename(ENSEMBL_ID = `Ensembl ID`,
                  HGNC_SYMBOL = `Gene Name`,
                  Barcode_Sequence = `Barcode sequence`,
                  Oligo_ID = DNA_ID,
                  Antigen = Description) %>%
    dplyr::mutate(Barcode_Sequence =
                      dplyr::coalesce(Barcode_Sequence, Barcode, Sequence),
                  Oligo_ID = dplyr::coalesce(Oligo_ID, `Format / Barcode`),
                  Clone = dplyr::coalesce(Clone, clone),
                  ENSEMBL_ID = dplyr::coalesce(ENSEMBL_ID, `Ensemble ID`),
                  HGNC_SYMBOL = dplyr::coalesce(HGNC_SYMBOL, `Gene name`),
                  Antigen = dplyr::coalesce(Antigen, Specificity)) %>%
    dplyr::select(-Barcode, -`Format / Barcode`, -clone,
                  -`Ensemble ID`, -`Gene name`, -Sequence, -Specificity) %>%
    dplyr::mutate(Antigen = gsub("mouse/human", "human/mouse", Antigen),
                  Reactivity =
                      AbNames:::.gsubNA("^anti-([Hh]uman(/mouse)?(/rat)?).*$",
                                  "\\1", Antigen),
                  Reactivity = tolower(Reactivity),
                  Antigen = gsub("^anti-([A-z\\/]+)\\s", "", Antigen),
                  # Two names start with Hu instead of "anti-human"
                  Antigen = gsub("^Hu\\s", "", Antigen),
                  Oligo_ID = substr(Oligo_ID, 2, nchar(Oligo_ID)),
                  Antigen = AbNames::replaceGreekSyms(Antigen, "sym2letter"))%>%
    dplyr::relocate(Antigen, Clone, ENSEMBL_ID, HGNC_SYMBOL, Oligo_ID,
                    TotalSeq_Cat, Barcode_Sequence, Reactivity)

# Fix an importing error
totalseq <- totalseq %>%
    dplyr::mutate(Antigen = ifelse(Antigen == "Fc?RI?", "FceRIA", Antigen))

# Some TotalSeq B Ensembl_IDs are duplicated barcode sequences, set to NA ----
totalseq <- totalseq %>%
    dplyr::mutate(ENSEMBL_ID = ifelse(grepl("^ENSG", ENSEMBL_ID),
                                      ENSEMBL_ID, NA))

# Check that there is only one Ensembl_ID per Antigen-Clone combo ----

totalseq <- totalseq %>%
    AbNames:::nPerGroup(group = c("Antigen", "Clone"), "ENSEMBL_ID")

if (max(totalseq$n_per_group) > 1){
    warning("More than one value of Ensembl_ID per Antigen/Clone combo")
}

totalseq <- totalseq %>% dplyr::select(-n_per_group)

# Fill in missing genes if Antigen, Clone, and Oligo_ID match ----
# Note that by adding "Reactivity", controls will be removed as reactivity was
# only defined for human
totalseq <- totalseq %>%
    AbNames::fillByGroup(group = c("Antigen", "Clone", "Reactivity"),
                         fill = c("ENSEMBL_ID", "HGNC_SYMBOL"))

# Fill in the missing gene symbols and reactivity if Ensembl ID is given
# (isotype controls do not have Ensembl IDs given)
totalseq <- totalseq %>%
    AbNames::fillByGroup(group = "ENSEMBL_ID",
                         fill = c("HGNC_SYMBOL", "Reactivity")) %>%
    dplyr::ungroup()

# Check that HGNC_SYMBOL in totalseq matches hgnc data set ----
data(hgnc)
hgnc <- hgnc %>%
    dplyr::select(HGNC_ID, HGNC_SYMBOL, ENSEMBL_ID) %>%
    dplyr::filter(! is.na(ENSEMBL_ID))

hgnc_to_ts <- hgnc$HGNC_SYMBOL[match(totalseq$ENSEMBL_ID, hgnc$ENSEMBL_ID)]
df <- tibble(hgnc = hgnc_to_ts, ts = totalseq$HGNC_SYMBOL)

# Values do not match when
# - totalseq is missing
# - totalseq gives extra information (TNFRSF8 (CD30), MME (CD10))
# - CD34 / CD334 - this is clearly a mistake in totalseq as Antigen = CD34

df %>% dplyr::filter(! hgnc == ts | is.na(ts))

# Replace symbols with HGNC values, add HGNC ids
totalseq <- totalseq %>%
    dplyr::select(-HGNC_SYMBOL) %>%
    dplyr::left_join(hgnc, by = "ENSEMBL_ID")

any(is.na(totalseq$HGNC_SYMBOL))

# Check the cases where multiple antigens are mapped to the same gene ----

x <- totalseq %>%
    dplyr::group_by(ENSEMBL_ID) %>%
    dplyr::mutate(n_antigen = n_distinct(Antigen)) %>%
    dplyr::filter(n_antigen > 1) %>%
    dplyr::select(ENSEMBL_ID, Antigen) %>%
    dplyr::arrange(ENSEMBL_ID) %>%
    unique()


# CD45 / CD45RA / CD45RO
# HLA-DR / HLA-DR, DP, DQ


# CD3 is CD3E on website but CD3D in totalseq, remove

totalseq <- totalseq %>%
    dplyr::filter(! Antigen == "CD3")



# Fill in missing ts_barcode Antigens using totalseq_cocktails -----
# ts_barcodes have NA for Antigen for isotype control.
# Fill in control name using totalseq_cocktails table
temp <- totalseq_cocktails %>%
    dplyr::select(Antigen, Clone) %>%
    unique() %>%
    group_by(Clone) %>%
    # If there is more than one antigen per clone, select the first
    # (checked manually, they are the same with different names)
    dplyr::slice_head(n = 1) %>%
    ungroup()

ts_barcodes <- ts_barcodes %>%
    dplyr::rows_patch(temp, unmatched = "ignore", by = c("Clone")) %>%
    dplyr::filter(! is.na(Antigen),
                  # When "Clone" is na it's Biotin
                  ! is.na(Clone))



# Create totalseq data set ----
totalseq <- as.data.frame(ts_barcodes)
usethis::use_data(totalseq, overwrite = TRUE, compress = "bzip2")


# Genes to check annotation: -----
# CD158 (KIR2DL1/S1/S3/S5)
# CD158b (KIR2DL2/L3, NKAT2)
# CD158e1 (KIR3DL1, NKB1)
# CD169
# CD3

# Exploration with inexact join (not necessary for this data set) -----

# Reactivity is more thoroughly described in ts_barcodes
# totalseq_cocktails often use longer antigen names

# Jaro Winkler distance is an edit distance with matches at the start given
# higher scores

# https://github.com/djvanderlaan/reclin

#library(reclin)
#
#p <- pair_blocking(totalseq_cocktails, ts_barcodes,
#                   blocking_var = c("Oligo_ID", "TotalSeq_Cat")) %>%
#    compare_pairs(by = c("Antigen", "Clone", "Oligo_ID", "Barcode_Sequence"),
#                  default_comparator = jaro_winkler(0.9), overwrite = TRUE) %>%
#    score_simsum(var = "simsum") %>%
#    select_greedy("simsum", var = "greedy", threshold = 0)
#
## Problem in problink_em - mprobs become 1 leading to div by zero
#
## For this data.set, it doesn't make a difference if greedy or n_to_m is used
#p <- data.frame( p[,c("x", "y")])




