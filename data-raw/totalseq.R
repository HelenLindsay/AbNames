# Notes ----

# TotalSeq_Barcodes is more complete than TotalSeq_Cocktails.
# TotalSeq_Cocktails just used to fill in missing isotype control names.

# Setup ----

library("readxl")
library("dplyr")
library("readr")
library("AbNames")
library("stringr")
library("janitor")
library("stringi")

#---------------------------------------------------------------------------
# Format TotalSeq barcode tables ----

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

# Filter to remove non-human gene identifiers
ts_barcodes <- ts_barcodes %>%
    splitUnnest(ab = "ENSEMBL_ID", split = ", ") %>%
    dplyr::mutate(ENSEMBL_ID = gsub("^NSG", "ENSG", ENSEMBL_ID), # missing E
                  ENSEMBL_ID = gsub("\\.[0-9]$", "", ENSEMBL_ID), # suffix
                  ENSEMBL_ID = ifelse(grepl("ENSG[0-9]+$", ENSEMBL_ID),
                                      ENSEMBL_ID, NA)) %>%
    dplyr::group_by(Cat_Number) %>%
    dplyr::mutate(ENSEMBL_ID = .toString(ENSEMBL_ID)) %>%
    unique()

# From this information, it appears that Antigen / Oligo / TotalSeq_Cat is
# enough to fill in Cat_Number and Clone
# Is Clone enough to fill in Antigen?

#---------------------------------------------------------------------------
# Format TotalSeq cocktail tables
#---------------------------------------------------------------------------
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

# Some TotalSeq B Ensembl_IDs are duplicated barcode sequences, set to NA ----
totalseq <- totalseq %>%
    dplyr::mutate(ENSEMBL_ID = ifelse(grepl("^ENSG", ENSEMBL_ID),
                                      ENSEMBL_ID, NA))

#---------------------------------------------------------------------------
# Fill barcodes using cocktails ----

# ts_barcodes have NA for Antigen for isotype control.
# Fill in Antigen name for controls using totalseq_cocktails table
temp <- totalseq %>%
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
                  ! is.na(Clone)) %>%
    dplyr::mutate(Cat_Number = as.character(Cat_Number),
                  ENSEMBL_ID = na_if(ENSEMBL_ID, "")) %>%
    dplyr::ungroup()

totalseq <- ts_barcodes

# Fix importing errors
imp_err <- c("TCR V&alpha;2" = "TCR Vα2", "Fc?RI?" = "FceRIA")

totalseq <- totalseq %>%
    dplyr::mutate(Antigen = stringr::str_replace_all(Antigen,
                                                     imp_err,
                                                     names(imp_err))) %>%
    # Select human or isotype control (Isotype controls have Reactivity NA)
    dplyr::filter(if_any(c(Reactivity, Cross_Reactivity), ~grepl("Human", .x)) |
                         (is.na(Reactivity) & is.na(ENSEMBL_ID)))

#---------------------------------------------------------------------------
# Checks
#---------------------------------------------------------------------------
# Check that there is only one Ensembl_ID per Antigen ----

# These were manually identified by anti_joining totalseq to hgnc,
# by looking at one antigen -> multiple genes or one gene -> multiple antigens,
# and checking the information on the biolegend website.n
# Note: "CDw" stands for "workshop", insufficiently characterised
fixes <- tibble::tribble(
    ~Antigen, ~Clone, ~ENSEMBL_ID,
    # Update Antigen with more information
    "CD209 (DC-SIGN)", "9E9A8", "ENSG00000090659",
    "CD209/CD299 (DC-SIGN/L-SIGN)", "14E3G7",
       "ENSG00000090659, ENSG00000104938",
    "CD85g", "17G10.2", "ENSG00000239961",
    "EGFR", "AY13", "ENSG00000146648",
    "GARP (LRRC32)", "7B11", "ENSG00000137507",
    "VEGFR-3 (FLT-4)", "9D9F9", "ENSG00000037280",
    "CD11a/CD18 (LFA-1)", "m24", "ENSG00000005844, ENSG00000160255",
    "CD207", "4C7", "ENSG00000116031",
    "TCR Vδ2", "B6", "ENSG00000211821",

    # www.beckman.ch, TCR Vα7.2 = TCRAV7S2, gene cards TCRAV7S2 = TRAV1-2
    # Gapin (2014) TRAV1-2 = TCR Vα7.2
    # From Biolegend: Vα7.2 joined with Jα33 is characteristic of MAIT cells
    "TCR Vα7.2", "3C10", "ENSG00000256553",

    # From totalseq website - clone SK1 recognises the alpha chain
    "CD8a", "SK1", "ENSG00000153563",

    # Add (NCAM) for consistency with other clone
    "CD56 (NCAM)", "QA17A16", "ENSG00000149294",

    # Overwrite incorrect ENSEMBL_ID
    # (discovered by multiple antigens to same gene id)
    # (verified via BioLegend / genenames websites)
    "CD6", "BL-CD6", "ENSG00000013725",
    "CD164", "67D2", "ENSG00000135535",
    "CD31", "WM59", "ENSG00000261371",
    "CD140a", "16A1", "ENSG00000134853",
    "CD24", "M1/69" ,"ENSG00000272398",
    "TCR Vγ9", "B3","ENSG00000211695",

    # Totalseq website: H1R2 reacts with common epitope of CD235a and CD235b
    "CD235ab", "HIR2", "ENSG00000170180, ENSG00000250361",

    # Totalseq website: CD98 is a heterodimer
    "CD98","MEM-108","ENSG00000103257, ENSG00000168003",

    # Totalseq website: 6D4 antibody reacts with common epitope of MICA/MICB
    "MICA/MICB", "6D4", "ENSG00000204516, ENSG00000204520",

    # Totalseq website: CD158 HP-MA4 rx with KIR2DL1, KIR2DS1, KIR2DS3, KIR2DS5
    # KIR2DS3 AND KIR2DS5 are on haplotype chromosomes
    "CD158", "HP-MA4",
        "ENSG00000125498, ENSG00000276387, ENSG00000277163, ENSG00000274739",

    # Totalseq website: CD158b DX27 reacts with common epitope of
    # KIR2DL2, KIR2DL3 and KIR2DS2
    "CD158b", "DX27", "ENSG00000274412, ENSG00000243772, ENSG00000278300",

    # Totalseq website: CD16 3G8 interacts with FcGRIIIa and FGγRIIIb receptors
    "CD16", "3G8", "ENSG00000203747, ENSG00000162747",

    # Totalseq website: CD66a/c/e binds epitope shared by CD66a, c and e
    "CD66a/c/e", "ASL-32", "ENSG00000079385, ENSG00000086548, ENSG00000105388",

    # Totalseq website: M1310G05 has higher affinity for
    # IgG1 and IgG3 than IgG2 and IgG4
    "IgG FC", "M1310G05", "ENSG00000211896, ENSG00000211897",

    "TCR Vα24-Jα18", "6B11", "ENSG00000211805, ENSG00000211871",

    # Totalseq website: Clone W6/32 reacts with beta2-microglobulin
    # (relabelling to beta2-microglobulin to avoid ambiguity with HLA gene ID)
    "HLA-A,B,C", "W6/32", "ENSG00000166710",

    # Assign CD3 to CD3E to match website
    "CD3", "UCHT1", "ENSG00000198851",
    "CD3", "SK7", "ENSG00000198851")


totalseq <- totalseq %>%
    dplyr::rows_update(fixes, by = "Clone", unmatched = "ignore")


# Now set ENSEMBL_ID to NA if there are multiple for the same antigen
totalseq <- totalseq %>%
    # Check for multiple ENSEMBL_IDs with the same Antigen
    AbNames:::nPerGroup(group = c("Antigen"), "ENSEMBL_ID") %>%
    # If there is more than one value of ENSEMBL_ID per Antigen set to zero
    # (18/7/22) - this is just CD3 mapped to CD3E and CD3D, fixed above
    dplyr::mutate(ENSEMBL_ID = ifelse(nENSEMBL_ID <= 1, ENSEMBL_ID, NA)) %>%
    dplyr::select(-nENSEMBL_ID)


# Fill in missing gene IDs if Antigen, Clone, and Oligo_ID match ----
totalseq <- totalseq %>%
    AbNames::fillByGroup(group = c("Antigen", "Clone"),
                         fill = c("ENSEMBL_ID"))

# Add HGNC_SYMBOL from hgnc data set ----
data(hgnc)
hgnc <- hgnc %>%
    dplyr::select(HGNC_ID, HGNC_SYMBOL, ENSEMBL_ID) %>%
    dplyr::filter(! is.na(ENSEMBL_ID)) %>%
    dplyr::mutate(SOURCE == "TotalSeq")
    unique()

totalseq <- left_join(totalseq, hgnc, by = c("ENSEMBL_ID"))


# Remove non-ascii characters -----

totalseq <- totalseq %>%
    dplyr::mutate(across(where(is.character),
                         ~stringi::stri_trans_general(.x,
                                    id="Any-Latin;Greek-Latin;Latin-ASCII")))


# Clean up column names -----





# Create totalseq data set ----
totalseq <- as.data.frame(totalseq)
usethis::use_data(totalseq, overwrite = TRUE, compress = "bzip2")


# Code for manually checking for inconsistencies ----

# # Antigens that have ENSEMBL ID in hgnc but do not share an alias
# aj <- semi_join(totalseq, hgnc, by = c("ENSEMBL_ID")) %>%
#     anti_join(hgnc, by = c("Antigen" = "value")) %>%
#     left_join(hgnc, by = "ENSEMBL_ID") %>%
#     dplyr::select(Antigen, HGNC_SYMBOL, value, Clone, ENSEMBL_ID) %>%
#     group_by(Antigen)
#
# # Check different antigen mapped to same gene
# sg <- totalseq %>%
#     group_by(ENSEMBL_ID) %>%
#     dplyr::filter(n_distinct(Antigen) > 1)
#
#
# # Select genes with null reactivity or Ensembl ID
# nrx <- totalseq %>%
#     dplyr::filter(is.na(Reactivity) | is.na(ENSEMBL_ID)) %>%
#     dplyr::select(-Barcode_Sequence, -Date_Released,
#                   -Oligo_ID, -TotalSeq_Cat) %>%
#     dplyr::arrange(Antigen, Clone) %>%
#     dplyr::group_by(Antigen, Clone)


# Notes ----

# Haven't overwritten but not sure
# TCR Va24 = TRAV10 Gapin "Check MAIT", https://www.beckman.ch
# (IGRa02 is only described sequence)
# https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:12103
# https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:12121

# from beckman TCR Vβ13.1 is also called TCRBV13S1
# TCRBV13S1 is HGNC alias for TRBV6-5 HGNC:12230 and TRBV13 HGNC:12188

# To do:
# CD45RA, CD45RB, CD45R0, ,TRAV24
#  HLA-A,B,C, HLA-A2
# HLA-DR L234 - does not cross react with DP and DQ - binds HLA-DRa
# HLA-DR, DP, DQ Tu39
# TRA-1-60-R != POXDL? TRA-1-81 != PODXL
# (from abcam) TRA-1-60 - a carbohydrate epitope associated with podocalyxin
# (from stemcell) TRA1-81 is a carbohydrate epitope associated with podocalyxin

