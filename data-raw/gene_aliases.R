# Notes -----

# Genes in NCBI not org_db are where ENTREZ_ID is mapped to a novel gene or
# scaffold.

# Both org.db and ncbi can map one ENTREZ_ID to several ENSEMBL_IDs

# NCBI contains Entrez IDs not in org_db and vice versa,
# usually when they differ it is for this reason

# Example of gene with ambiguous aliases -
# CCL4L1/2 both have alias AT744.2

# NCBI genes has more cases where one HGNC symbol is mapped to several
# ENSEMBL IDs

# Prefer to use two matches (e.g. symbol + ensembl) if possible

# Not all HGNC Ensembl IDs are in the biomaRt table
# (some are not protein-coding)

# ENSG00000203668: Entrez gene symbol CHML, HGNC symbol OPN3 -
# both are valid HGNC symbols

# Within hgnc, one Uniprot ID can map to multiple HGNC IDs, e.g. Q9GZY0

# Setup -----

library(tidyverse)
library(readxl)

source("ncbi.R")
source("org_db.R")


# Add Ensembl_ID to HGNC where Ensembl matches but HGNC doesn't  -----



# Merge HGNC, NCBI and org.db ---------------------------------------

prep_merge <- function(df, hgnc){
    df %>%
        dplyr::select(ENSEMBL_ID, ENTREZ_ID, HGNC_SYMBOL) %>%
        dplyr::semi_join(hgnc) %>%
        # Want Ensembl, Entrez and Symbol for merging
        na.omit() %>%
        unique()
}

common_or_missing <- function(x, y){
    sj <- dplyr::semi_join(x, y)
    aj <- dplyr::anti_join(x, y)
    # Are the entries in the anti-join inconsistent or missing?

    # Find entries in y where any column value is in the anti_join
    uj <- union_join(y, aj)

    # Find entries in the anti-join where any entry is in y (using above)
    to_rm <- union_join(aj, uj)

    aj <- anti_join(aj, to_rm)
    return(dplyr::bind_rows(sj, aj))
}


bm_x <- prep_merge(bm, hgnc)
ncbi_x <- prep_merge(ncbi_genes, hgnc)
org_db_x <- prep_merge(org_db, hgnc)

bm_ncbi <- common_or_missing(bm_x, ncbi_x)
bm_orgdb <- common_or_missing(bm_x, org_db_x)
ncbi_bm <- common_or_missing(ncbi_x, bm_x)
ncbi_orgdb <- common_or_missing(ncbi_x, org_db_x)
orgdb_bm <- common_or_missing(org_db_x, bm_x)
orgdb_ncbi <- common_or_missing(org_db_x, ncbi_x)

consistent <- Reduce(dplyr::full_join,
                     list(bm_ncbi, bm_orgdb, ncbi_bm, ncbi_orgdb,
                                  orgdb_bm, orgdb_ncbi))

hgnc <- hgnc %>% dplyr::left_join(consistent)


# Add novel aliases ------------------------------------------------

# Add novel aliases from biomaRt and NCBI
# (Note - have previously checked that HGNC_ID / ENSEMBL_ID combinations
# are correct)
hgnc <- hgnc %>%
    dplyr::bind_rows(bm_novel, ncbi_novel, org_db_novel) %>%

    # Fill Entrez IDs from new aliases
    dplyr::group_by(HGNC_ID) %>%
    tidyr::fill(ENTREZ_ID, HGNC_NAME, UNIPROT_IDS, .direction = "updown") %>%

    # Fill HGNC IDs (missing from org_db)
    dplyr::group_by(ENSEMBL_ID) %>%
    tidyr::fill(HGNC_ID, HGNC_NAME, UNIPROT_IDS, .direction = "updown") %>%

    # New values may come from more than one source, aggregate
    dplyr::group_by(HGNC_ID, ENSEMBL_ID, UNIPROT_IDS, HGNC_SYMBOL, ENTREZ_ID,
                    symbol_type, value) %>%
    dplyr::summarise(SOURCE = toString(SOURCE), .groups = "keep") %>%

    # Check if values are unambiguous
    # (newly added aliases may create new ambiguities, only 3 are ambiguous)
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_ID)) %>%
    dplyr::filter(n_genes == 1 | SOURCE == "HGNC") %>%
    dplyr::select(-n_genes)


#> table(is.na(hgnc$ENTREZ_ID)) # Note - not the number of genes because aliases
# 122569  36302

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
