# Cell marker database ---------------------------------------------

# Note: Isoform name may be mapped to gene name, e.g. CD45RO -> PTPRC
# Some antigens have no annotation information, e.g. CD45.1

# http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.jsp

gsubCellmarker <- function(x){
    p1_matches <- stringr::str_extract_all(x, "\\[([^\\[]+)\\]")
    p1_sub <- gsub(", ", "\\|", unlist(p1_matches))
    p1_sub <- gsub("\\[|\\]", "", p1_sub)
    p1_sub <- relist(p1_sub, p1_matches)
    m <- gregexpr("\\[([^\\[]+)\\]", x)
    regmatches(x, m) <- p1_sub
    return(x)
}

# Cellmarker
# Sometimes family name is listed under "ENTREZ_IDS"
#

# To do: make sure that cellmarker cellMarker col mapped to same geneSymbol/ID
# column is a correct, unambiguous alias
# Is ENTREZ symbol the HGNC symbol? - No, sometimes it's previous symbol
# or not found

cellmarker_loc <- paste0("http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/",
                         "Human_cell_markers.txt")
cellmarker_fname <- "~/Analyses/CITEseq_curation/data/CellMarker_human.txt"
download.file(cellmarker_loc, destfile = cellmarker_fname)

cellmarker <- readr::read_delim(cellmarker_fname) %>%
    dplyr::select(cellName, cellMarker, geneSymbol,
                  geneID, proteinName, proteinID) %>%
    dplyr::filter(if_all(c(geneID, proteinID), ~! is.na(.x))) %>%
    unique() %>%
    # Separate protein complexes with | instead of [a, b, c]
    dplyr::mutate(across(c(geneSymbol, geneID, proteinName, proteinID),
                         ~ gsubCellmarker(.x))) %>%
    dplyr::mutate(across(c(cellMarker, geneSymbol, geneID,
                           proteinName, proteinID), ~strsplit(.x, ", "))) %>%

    # If there is not 1 marker per gene, something is wrong
    dplyr::filter(lengths(cellMarker) == lengths(geneSymbol) &
                      lengths(cellMarker) == lengths(proteinID) &
                      lengths(geneSymbol) == lengths(geneID) &
                      lengths(proteinID) == lengths(proteinName)) %>%
    tidyr::unnest(cols = c(cellMarker, geneSymbol, geneID,
                           proteinName, proteinID)) %>%
    dplyr::rename(Antigen = cellMarker,
                  ENTREZ_IDS = geneID,
                  ENTREZ_SYMBOL = geneSymbol,
                  UNIPROT_IDS = proteinID) %>%
    dplyr::mutate(SOURCE = "CELLMARKER") %>%
    dplyr::mutate(across(c(ENTREZ_SYMBOL, ENTREZ_IDS,
                           proteinName, UNIPROT_IDS), ~na_if(., "NA"))) %>%
    dplyr::mutate(across(everything(), ~stringr::str_squish(.x)))


# Exploration ----
# cellmarker %>% dplyr::select(-cellName) %>% unique()
# 11,231 Antigens, 1,250 don't match (exactly) via symbol
# (excluding families, etc)
# Why?
# Because it's long form
# e.g. "Calcitonin and Related Receptor Antagonists"
# Because it's not very specific:
# Stage Specific Embryonic Antigens
# Because it's an RNA - C1orf61 (!)
# Because it matches a previous symbol
# Because the antigen is the official symbol but there is no gene info
# Because it's a pseudogene! EGFEM1P

x <- cellmarker %>%
    dplyr::select(-cellName) %>%
    unique() %>%
    dplyr::filter(! ENTREZ_SYMBOL %in% hgnc$HGNC_SYMBOL, !
                      grepl("\\|", ENTREZ_SYMBOL),
                  ! grepl("family", ENTREZ_SYMBOL))


# -----





protein_complexes <- cellmarker %>%
    dplyr::filter(grepl("\\|", ENTREZ_SYMBOL)) %>%
    dplyr::select(Antigen, ENTREZ_SYMBOL, ENTREZ_IDS,
                  proteinName, UNIPROT_IDS) %>%
    unique()

# Using org.db, add ENSEMBL identifiers (one ENSEMBL may map to several ENTREZ?)
# As this includes "ALIAS" column, rows may be added
protein_complexes <- protein_complexes %>%
    dplyr::mutate(across(c(ENTREZ_SYMBOL, ENTREZ_IDS,
                           proteinName, UNIPROT_IDS), ~strsplit(.x, "\\|"))) %>%
    tidyr::unnest(cols = c(ENTREZ_SYMBOL, ENTREZ_IDS,
                           proteinName, UNIPROT_IDS)) %>%

    dplyr::left_join(org_db %>% dplyr::select(-value) %>% unique(),
                     by = c(ENTREZ_IDS = "ENTREZ_ID")) %>%

    dplyr::group_by(Antigen) %>%



    # When a single gene is mapped to several Antigens, check if cellName annotation
    # is also the same

    # Split into individual genes, check the symbols are official symbols,
    # add HGNC and ENSEMBL IDs.





# Row Epithelial cell starting Adhesion molecules -
# Number of cell markers is not equal to number of gene (groups),
# is there an extra , Adhesion molecules, LFA1, Adhesion molecules LFA2?

# Cell surface protein atlas ---------------------------------------

# http://wlab.ethz.ch/cspa/#abstract

# Notes:
# Putative Ig-like domain-containing protein / 374383
# ODZ3 / 55714 = TENM3 checked manually via gene cards

# Most entries can be matched to HGNC via UNIPROT_ID / HGNC_SYMBOL
# Those that can't include HLA and T-cell receptors

# Does not appear to contain protein complex names


# Table of validated surfaceome proteins
#cspa_loc <- "http://wlab.ethz.ch/cspa/data/S2_File.xlsx"
#cspa_fname <- "~/Analyses/CITEseq_curation/data/cspa.xlsx"
#download.file(cspa_loc, destfile = cspa_fname)

# Sheet A has human proteins and entrez IDs
cspa <- #readxl::read_xlsx(cspa_fname) %>%
    readxl::read_xlsx("~/Analyses/CITEseq_curation/data/cspa.xlsx") %>%
    dplyr::select(UP_Protein_name, CD, ENTREZ_gene_ID,
                  `ENTREZ gene symbol`, ID_link) %>%
    dplyr::rename(UNIPROT_NAME = UP_Protein_name,
                  UNIPROT_ID = ID_link,
                  ENTREZ_ID = ENTREZ_gene_ID,
                  ENTREZ_SYMBOL = `ENTREZ gene symbol`,
                  Antigen = CD) %>%
    dplyr::mutate(SOURCE = "CSPA",
                  Antigen = na_if(Antigen, "no"),
                  across(c(ENTREZ_SYMBOL, ENTREZ_ID), ~na_if(.x, "0")),
                  ENTREZ_ID = as.character(ENTREZ_ID)) %>%
    # Join by exact match to "value" in HGNC
    # (either official symbol or alias)
    dplyr::left_join(hgnc %>%
                         dplyr::select(HGNC_SYMBOL, HGNC_ID, value) %>%
                         unique(),
                     by = c("ENTREZ_SYMBOL" = "value"))

# Patch via UNIPROT_ID, if it is unique within HGNC
up_ids <- cspa %>%
    dplyr::filter(is.na(HGNC_ID)) %>%
    dplyr::pull(UNIPROT_ID)

hgnc_patch <- hgnc %>%
    dplyr::filter(UNIPROT_ID %in% up_ids) %>%
    dplyr::group_by(UNIPROT_ID) %>%
    dplyr::filter(n_distinct(HGNC_ID) == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(HGNC_SYMBOL, HGNC_ID, UNIPROT_ID) %>%
    unique()

cspa <- cspa %>%
    dplyr::rows_patch(hgnc_patch, by = "UNIPROT_ID") %>%
    # 5 unmatched have only UNIPROT IDs, some are obsolete
    dplyr::filter(! is.na(HGNC_ID))

# Reformat cspa for joining to hgnc
cspa <- cspa %>%
    dplyr::mutate("CSPA" = TRUE,
                  ENTREZ_SYMBOL =
                      AbNames:::.noDups(ENTREZ_SYMBOL, HGNC_SYMBOL)) %>%
    dplyr::rename(ALIAS = ENTREZ_SYMBOL) %>%
    tidyr::pivot_longer(cols = c("UNIPROT_NAME", "Antigen", "ENTREZ_SYMBOL"),
                        values_to = "value", names_to = "symbol_type") %>%
    dplyr::filter(! is.na(value))

# Check for inconsistencies -----
# TO DO: HAVEN'T JOINED IN ENTREZ YET!

cspa %>% dplyr::anti_join(hgnc,
                          by = c("HGNC_ID", "ENTREZ_ID", "UNIPROT_ID"))


# Add a column to hgnc indicating if gene is in CSPA
hgnc <- hgnc %>%
    dplyr::left_join(cspa %>% dplyr::select(HGNC_ID, CSPA))

cspa_novel <- cspa %>%
    dplyr::anti_join(hgnc, by = c("HGNC_ID", "HGNC_SYMBOL",
                                  "UNIPROT_ID", "value"))
