library("tidyverse")

# Download from antibody registry requires a login

ar_all <- "inst/extdata/antibody_registry_2023-05-15.csv"
ar_human <- "inst/extdata/antibody_registry_human_monoclonal_2023-05-15.csv"
ar <- read_delim(ar_all) %>%
    dplyr::filter(grepl("Anti-Human", ab_name),
                  grepl("Monoclonal", ab_name),
                  ! grepl("Fluorescein|PE Conju", ab_name),
                  grepl("[Cc]lone", ab_name))

readr::write_csv(ar, file = ar_human)

ar <- ar %>%
    # Fix formatting errors
    dplyr::mutate(ab_name = gsub("MouseAnti", "Mouse Anti", ab_name),
                  ab_name = gsub("Clone [Cc]lone", "Clone", ab_name),
                  ab_name = gsub("AntibodyMouse.*", "Antibody", ab_name)) %>%
    tidyr::separate(ab_name, into = c("Species", "ab_name"),
                    sep = " Anti-Human", fill = "left") %>%
    dplyr::mutate(across(where(is_character), stringr::str_squish)) %>%
    tidyr::separate(ab_name, into = c("ab_name", "Clone"),
                    sep = "[,]? ?[Cc]lone[s]? ",
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
