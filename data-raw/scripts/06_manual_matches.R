mm <- read_csv("inst/extdata/rols_ontology.csv") %>%
    dplyr::mutate(across(where(is_character), stringr::str_squish))
mm_ids <- mm %>% dplyr::select(Antigen, Clone, HGNC_ID, HGNC_SYMBOL, ALT_ID)

# Add ALT_ID to gene aliases data set -----

# (for regenerating)
gene_aliases <- gene_aliases %>%
    dplyr::filter(! SOURCE == "MANUAL_LOOKUP") %>%
    dplyr::select(-any_of("ALT_ID"))

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

# Create gene_aliases data set
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
