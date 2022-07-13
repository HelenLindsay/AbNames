library(tidyverse)

# Download from antibody registry requires a login

ar_all <- "~/Analyses/CITEseq_curation/data/antibody_registry.csv"
ar_human <- paste0("~/Analyses/CITEseq_curation/data/",
                   "antibody_registry_human_monoclonal.csv")
ar <- read_delim(ar_all) %>%
    dplyr::filter(grepl("Anti-Human", ab_name),
                  grepl("Monoclonal", ab_name),
                  ! grepl("Fluorescein|PE Conju", ab_name),
                  grepl("[Cc]lone", ab_name))

readr::write_delim(ar, file = ar_human)

ar <- ar %>%
    # Fix formatting errors
    dplyr::mutate(ab_name = gsub("MouseAnti", "Mouse Anti", ab_name),
                  ab_name = gsub("Clone [Cc]lone", "Clone", ab_name),
                  ab_name = gsub("AntibodyMouse.*", "Antibody", ab_name)) %>%
    tidyr::separate(ab_name, into = c("Species", "ab_name"),
                    sep = " Anti-Human", fill = "left") %>%
    dplyr::mutate(across(where(is_character), stringr::str_squish)) %>%
    tidyr::separate(ab_name, into = c("ab_name", "Clone"), sep = "[,]? ?[Cc]lone[s]? ",
                    remove = FALSE, extra = "merge") %>%
    dplyr::select(-defining_citation) %>%
    dplyr::mutate(Antigen = gsub(" Monoclonal.*", "", ab_name)) %>%
    dplyr::filter(! ab_name == Antigen,
                  ! grepl("[Cc]lone", Clone)) %>%
    dplyr::rename(RRID = id,
                  Citation = proper_citation,
                  Cat_Number = catalog_num,
                  Vendor = vendor) %>%
    dplyr::select(-ab_name)
