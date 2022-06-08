# Notes -----

# Genes in NCBI not org_db are where ENTREZ_ID is mapped to a novel gene or
# scaffold.

# Both org.db and ncbi can map one ENTREZ_ID to several ENSEMBL_IDs

# NCBI contains Entrez IDs not in org_db

# To do ----

# Check gene biotype of ensembl genes

# Check ensembl mapping location matches map_location (only one should?)

# check whether each ensembl id in hgnc maps to only one hgnc symbol / id

# check whether unofficial symbols map to a different official gene
# ENSG00000203668: Entrez gene symbol CHML, HGNC symbol OPN3 -
# both are valid HGNC symbols

# Single NCBI maps to multiple ENSEMBL, are both locations valid?
# Do they produce the same protein?  Does Protein Ontology differentiate?

# Setup -----

library(tidyverse)
library(readxl)

source("ncbi.R")
source("org_db.R")


# Merge HGNC, NCBI and org.db --------------------------------------

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
    dplyr::group_by("HGNC_ID", "value") %>%
    dplyr::mutate(value = toString(SOURCE)) %>%
    unique() %>%

    # Check if values are unambiguous
    # (newly added aliases may create new ambiguities)
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_ID)) %>%
    dplyr::filter(n_genes == 1 | SOURCE == "HGNC") %>%
    dplyr::select(-n_genes)









# Check if ncbi_genes and org_db agree on mapping between
# HGNC, ENTREZ and ENSEMBL -----

ncbi_id_map <- ncbi_genes %>%
    dplyr::select(NCBI_ID, ENSEMBL_ID) %>%
    unique()

org_db_id_map <-  org_db %>%
    dplyr::select(ENTREZ_ID, ENSEMBL_ID) %>%
    unique()

ncbi_diff <- anti_join(ncbi_id_map, org_db_id_map)
org_db_diff <- anti_join(org_db_id_map, ncbi_id_map)



# Joining Entrez and Ensembl -----

# Usually one ENSEMBL to one HGNC (occasional exceptions), e.g. TCBE
# TCBE appears to be multiple copies that produce the same protein
# ENTREZ can be several to one HGNC

# Within org_db, one ENTREZ_ID can map to several ENSEMBL_IDs
# e.g. APOBEC3A_B / APOBEC3A BiomaRt only maps to one of the aliases



# Sometimes biomaRt ENTREZ_IDs are NA.
# Fill these if the ENSEMBL_ID and HGNC_SYMBOL agree

org_db_patch <- org_db %>%
    dplyr::select(ENSEMBL_ID, HGNC_SYMBOL, ENTREZ_ID) %>%
    dplyr::filter(! is.na(ENSEMBL_ID)) %>%
    unique()

bm <- bm %>%
    dplyr::rows_patch(org_db_patch,
                      by = c("ENSEMBL_ID", "HGNC_SYMBOL"),
                      unmatched = "ignore")






union_join <- function(df, rows){
    x <- unlist(df[rows, ], use.names = FALSE)
    df %>% dplyr::filter(if_any(.cols = everything(), ~.x %in% x))
}




# There are cases where Biomart and org_db map the same
# ENSEMBL_ID / HGNC_SYMBOL to different ENTREZ_IDs,
# e.g. ENSG00000143702/CEP170 to 645455 (biomaRt) or 9859 (org_db)
# In this case, biomaRt mapping is to a pseudogene according to NCBI


# Just keep rows from org_db and Biomart where Entrez / Ensembl agree

bm_minus_entrez <- dplyr::anti_join(bm, org_db,
                                    by = c("ENSEMBL_ID", "ENTREZ_ID"))

ens_entrez <- dplyr::semi_join(bm, org_db,
                               by = c("ENSEMBL_ID", "ENTREZ_ID"))


# Two problems: one is mapping the IDs, the other is collecting the aliases




# Check for disagreements between org_db and NCBI ----
ncbi_not_org_db <- ncbi_genes %>%
    dplyr::select(ENTREZ_ID, ENSEMBL_ID) %>%
    unique() %>%
    dplyr::anti_join(org_db %>%
                         dplyr::select(ENTREZ_ID, ENSEMBL_ID) %>%
                         unique())
































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
# To do: make sure that cellmarker cellMarker col mapped to same geneSymbol/ID
# column is a correct, unambiguous alias

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
                  UNIPROT_IDS = proteinID)


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
    dplyr::left_join(org_db %>% dplyr::select(-ALIAS) %>% unique(),
                     by = c(ENTREZ_IDS = "ENTREZ_ID"))



# When a single gene is mapped to several Antigens, check if cellName annotation
# is also the same



# Split into individual genes, check the symbols are official symbols,
# add HGNC and ENSEMBL IDs.





# Row Epithelial cell starting Adhesion molecules -
# Number of cell markers is not equal to number of gene (groups),
# is there an extra , Adhesion molecules, LFA1, Adhesion molecules LFA2?


#dplyr::rename(ENTREZ_ID = geneID,
#             UNIPROT_ID = proteinID) %>%


# Cell surface protein atlas ---------------------------------------

# http://wlab.ethz.ch/cspa/#abstract


# Table of validated surfaceome proteins
cspa_loc <- "http://wlab.ethz.ch/cspa/data/S2_File.xlsx"
cspa_fname <- "~/Analyses/CITEseq_curation/data/cspa.xlsx"
download.file(cspa_loc, destfile = cspa_fname)

# Sheet A has human proteins and entrez IDs
cspa <- readxl::read_xlsx(cspa_fname) %>%
    dplyr::select(UP_Protein_name, UP_entry_name, CD, ENTREZ_gene_ID,
                  `ENTREZ gene symbol`) %>%
    dplyr::rename(UNIPROT_NAME = UP_Protein_name,
                  ENTREZ_ID = ENTREZ_gene_ID,
                  ENTREZ_SYMBOL = `ENTREZ gene symbol`,
                  Antigen = CD)



