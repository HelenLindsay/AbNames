library("tidyverse")
library("AbNames")

data(totalseq)
data(gene_aliases)

# Add ALT_ID to totalseq data set -----

original_nrow <- nrow(totalseq)

mm_fname <- system.file("extdata", "rols_ontology.csv", package = "AbNames")
mm <- read_delim(mm_fname) %>%
    dplyr::mutate(across(where(is_character), stringr::str_squish))
mm_ids <- mm %>% dplyr::select(Antigen, Clone, HGNC_ID, HGNC_SYMBOL, ALT_ID)

# We do not want to join by HGNC_ID as e.g. TRA-1-60-R matches PODXL but is not
# in the totalseq table
totalseq <- totalseq %>%
    left_join_any(mm_ids, cols = c("Antigen", "Clone"), shared = "update") %>%
    dplyr::mutate(ALT_ID = dplyr::coalesce(ALT_ID, HGNC_ID))

new_nrow <- nrow(totalseq)

# Check that totalseq did not gain rows during merge
original_nrow == new_nrow


totalseq <- as.data.frame(totalseq)
usethis::use_data(totalseq, overwrite = TRUE, compress = "bzip2")

# Manual check that all totalseq genes have a match
# totalseq %>%
# dplyr::filter(is.na(HGNC_ID), grepl("[Hh]uman", Reactivity)) %>%
# select(Antigen, Clone) %>%
# filter(! (Antigen %in% mm$Antigen |
# Clone %in% mm$Clone | Antigen %in% gene_aliases$value)) %>%
# unique() %>% arrange(Antigen) %>% data.frame()

# Add ALT_ID to gene aliases data set -----

mm_ids <- mm_ids %>%
    dplyr::select(-Clone) %>%
    unique() %>%
    dplyr::rename(value = Antigen)  %>%
    dplyr::mutate(SOURCE = "MANUAL_LOOKUP")

gene_aliases <- gene_aliases %>%
    dplyr::full_join(mm_ids, by = c("value", "HGNC_ID", "HGNC_SYMBOL")) %>%
    dplyr::mutate(ALT_ID = dplyr::coalesce(ALT_ID, HGNC_ID),
                  SOURCE = dplyr::coalesce(`SOURCE.y`, `SOURCE.x`)) %>%
    dplyr::select(-`SOURCE.x`, -`SOURCE.y`)

# Recreate gene_aliases data set
gene_aliases <- as.data.frame(gene_aliases)
usethis::use_data(gene_aliases, overwrite = TRUE, compress = "bzip2")


# Antibody notes -----

# CD158f/ KIR2DL5
# https://doi.org/10.3389/fimmu.2012.00289
# KIR2DL5 is highly polymorphic and exhibits copy number variation
# KIR2DL5A (telomeric) and KIR2DL5B (centromeric, usually inactivated)
# KIR2DL5 is represented in the Immuno Polymorphism Database by 15 KIR2DL5A
# and 25 KIR2DL5B alleles.  One single polymorphism distinguishes A and B

#IFN_gamma <- paste("https://proconsortium.org/app/entry/PR:000001361/",
#                   "IFN-gamma receptor 1 and IFN-gamma-R-alpha are synonyms")

#HLA_DR <- paste("Grep in protein ontology for HLA-DR gives",
#                "wrong results PR:P04233, PR:000001822")

#HLA_ABC <- paste("Biolegend website: Clone W6/32 recognizes ",
#                 "human beta2-microglobulin")

#CD66 <- paste0("https://www.bdbiosciences.com/en-ch/products/reagents/",
#               "flow-cytometry-reagents/research-reagents/",
#               "single-color-antibodies-ruo/pe-mouse-anti-human-cd66.551480")

#CD158 <- paste0("From BioLegend website: ",
#                "mAb HP-MA4 reacts with KIR2DL1 (CD158a), ",
#                "KIR2DS1 (CD158h), KIR2DS3, and KIR2DS5 (CD158g).")
#IgG <- paste0("From BioLegend website: clone M1310G05 has stronger affinity ",
#              "for IgG1 and IgG3 than for IgG2 and IgG4")

#TRAV7 <- paste0("Vα7.2 TCR with Jα33 forms invariant T cell receptor")

#CLA <- paste0("From BioLegend website: ",
#              "CLA is a scarbohydrate epitope of sialic acid ",
#              "and fucose-modified P-selectin glycoprotein ligand-1 (PSGL-1) ",
#              "Ligand for E-selectin, P-selectin, and L-selectin.")
