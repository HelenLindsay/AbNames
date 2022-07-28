data(gene_aliases)

# Cell marker database ---------------------------------------------

# Note: Isoform name may be mapped to gene name, e.g. CD45RO -> PTPRC
# Some antigens have no annotation information, e.g. CD45.1
# Table contains annotation errors, e.g. row shift error in UNIPROT_IDs

# Row Epithelial cell starting Adhesion molecules -
# Number of cell markers is not equal to number of gene (groups),
# is there an extra , Adhesion molecules, LFA1, Adhesion molecules LFA2?

# http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.jsp

# For patching ENTREZ_SYMBOL using Antigen
ga_patch <- gene_aliases %>%
    dplyr::select(value, UNIPROT_ID, ENTREZ_ID) %>%
    dplyr::rename(Antigen = value) %>%
    unique()

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

# To do: make sure that cellmarker cellMarker col mapped to same geneSymbol/ID
# column is a correct, unambiguous alias
# Is ENTREZ symbol the HGNC symbol? - No, sometimes it's previous symbol
# or not found

#cellmarker_loc <- paste0("http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/",
#                         "Human_cell_markers.txt")
cellmarker_fname <- "~/Analyses/CITEseq_curation/data/CellMarker_human.txt"
#download.file(cellmarker_loc, destfile = cellmarker_fname)


cellmarker <- readr::read_delim(cellmarker_fname) %>%
    dplyr::select(cellName, cellMarker, geneSymbol,
                  geneID, proteinName, proteinID) %>%
    dplyr::filter(if_all(c(geneID, proteinID), ~! is.na(.x)),
                  # "x family" isn't useful as an ID
                  ! grepl("family", geneID)) %>%
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
    # One row per gene
    tidyr::unnest(cols = c(cellMarker, geneSymbol, geneID,
                           proteinName, proteinID)) %>%
    dplyr::rename(Antigen = cellMarker,
                  ENTREZ_ID = geneID,
                  ENTREZ_SYMBOL = geneSymbol,
                  UNIPROT_ID = proteinID) %>%
    dplyr::mutate(SOURCE = "CELLMARKER") %>%
    dplyr::mutate(across(c(ENTREZ_SYMBOL, ENTREZ_ID,
                           proteinName, UNIPROT_ID), ~na_if(., "NA"))) %>%
    # Patch missing ENTREZ or UNIPROT ids
    dplyr::coalesce(ga_patch)




    # Filter again, some cellNames include an NA in a list of markers
    # (some aren't protein-coding, some would have an ENTREZ ID)
    dplyr::filter(if_all(c(ENTREZ_ID, UNIPROT_ID), ~! is.na(.x))) %>%


    dplyr::mutate(across(everything(), ~stringr::str_squish(.x))) %>%
    # Split the protein complexes
    dplyr::mutate(across(c(ENTREZ_SYMBOL, ENTREZ_ID, proteinName, UNIPROT_ID),
                         ~strsplit(.x, "\\|"))) %>%
                               tidyr::unnest(cols = c(ENTREZ_SYMBOL,
                                                      ENTREZ_ID,
                                                      proteinName,
                                                      UNIPROT_ID)) %>%
    # Not useful if antigen already exists
    dplyr::filter(! Antigen == ENTREZ_SYMBOL,
                  ! Antigen %in% gene_aliases$value)

# There are 13 markers that don't match HGNC_SYMBOL / ENTREZ_ID combo
# Removing, but probably safe to keep them
cellmarker <- cellmarker %>%
    dplyr::semi_join(gene_aliases,
                     by = c("ENTREZ_SYMBOL" = "HGNC_SYMBOL", "ENTREZ_ID")) %>%
    # Patch the uniprot IDs
    dplyr::rename(HGNC_SYMBOL = ENTREZ_SYMBOL) %>%
    dplyr::rows_update(hgnc %>%
                           dplyr::select(HGNC_SYMBOL, ENTREZ_ID, UNIPROT_ID) %>%
                           unique(),
                       unmatched = "ignore") %>%

    # Re-aggregate the protein complexes
    unique() %>%
    dplyr::group_by(Antigen, cellName) %>%
    dplyr::mutate(across(c(HGNC_SYMBOL, ENTREZ_ID, proteinName, UNIPROT_ID),
                         ~ paste(.x, collapse = "|"))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-cellName) %>%
    unique() %>%

    # Protein names will be added as an alias if they don't already exist
    dplyr::mutate(proteinName = ifelse(proteinName == HGNC_SYMBOL,
                                       NA, proteinName),
                  proteinName = ifelse(proteinName %in% gene_aliases$value,
                                       NA, proteinName)) %>%

    # Pivot longer for merging with aliases
    dplyr::rename(CELLMARKER_ANTIGEN = Antigen,
                  CELLMARKER_PROTEIN = proteinName) %>%
    tidyr::pivot_longer(cols = c("CELLMARKER_ANTIGEN", "CELLMARKER_PROTEIN"),
                        names_to = "symbol_type", values_to = "value") %>%
    # Aggregated protein names probably aren't useful
    dplyr::filter(! is.na(value), ! grepl("\\|", value))

# To do: check for same antigen to multiple genes (e.g. CD45RO)



# Exploration ----
# cellmarker %>% dplyr::select(-cellName) %>% unique()
# 11,231 Antigens, 1,250 don't match (exactly) via symbol
# (excluding families, etc)
# Why?
# Because it's long form
# e.g. "Calcitonin and Related Receptor Antagonists"
# Because it's not very specific:
# Stage Specific Embryonic Antigens
# Because it's an RNA - C1orf61
# Because it matches a previous symbol
# Because the antigen is the official symbol but there is no gene info
# Because it's a pseudogene! EGFEM1P

# Demonstration of row shift error in UNIPROT_IDs
#missing <- anti_join(cellmarker %>% dplyr::mutate(row_n = row_number()),
#                     gene_aliases, by = c("ENTREZ_ID", "UNIPROT_ID"))
#x <- missing$row_n - seq_along(missing$row_n)
#rle(x)$lengths

# Notes:
# Possible to get e.g. CD32, CD16, CD3 from other source?
# Haven't checked if individual genes are filtered from protein complexes

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
    dplyr::filter(! is.na(ENTREZ_SYMBOL)) %>%
    # Join by exact match to "value" in gene_aliases table
    # (either official symbol or alias)
    dplyr::inner_join(gene_aliases %>%
                          dplyr::select(HGNC_SYMBOL, HGNC_ID,
                                        ENTREZ_ID, value) %>%
                          unique(),
                      by = c("ENTREZ_SYMBOL" = "value", "ENTREZ_ID")) %>%
    # Know from join that ENTREZ_SYMBOL must match an existing symbol or alias
    dplyr::mutate(Antigen = ifelse(Antigen %in%
                                       gene_aliases$value, NA, Antigen)) %>%
    dplyr::filter(! is.na(Antigen)) %>%
    dplyr::select(HGNC_ID, HGNC_SYMBOL, ENTREZ_ID,
                  UNIPROT_ID, Antigen, SOURCE) %>%
    dplyr::rename(value = Antigen)


# Check for inconsistencies -----
# Two entries from CSPA have Uniprot ID, NA in HGNC

cspa %>% dplyr::anti_join(gene_aliases,
                          by = c("HGNC_ID", "ENTREZ_ID", "UNIPROT_ID"))



# Add a column to hgnc indicating if gene is in CSPA
gene_aliases <- gene_aliases %>%
    dplyr::left_join(cspa %>% dplyr::select(HGNC_ID, CSPA))

cspa_novel <- cspa %>%
    dplyr::anti_join(gene_aliases,
                     by = c("HGNC_ID", "HGNC_SYMBOL", "UNIPROT_ID", "value"))
