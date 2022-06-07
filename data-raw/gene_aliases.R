# To do ----

# Check gene biotype of ensembl genes
# Remove if same symbol mapped to two different ensembl genes?
# e.g. one NCBI SYMBOL to many HGNC symbols or ENSEMBL genes
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

# Merge HGNC and NCBI Entrez gene -----

# Only want genes for which there is a HGNC / ENSEMBL combination in HGNC
# (Note that not all HGNC IDs are in hgnc data set, e.g. non-protein-coding)

# There are some differences in which ENSEMBL_ID is mapped to which HGNC_ID,
# e.g. in HGNC HGNC:4883 -> ENSG00000000971 (official gene)
#      in NCBI HGNC:4883 -> ENSG00000289697 (novel gene)

ncbi_f <- ncbi_genes %>%
    dplyr::semi_join(hgnc %>% dplyr::select(HGNC_ID, ENSEMBL_ID))

# NCBI and HGNC can differ in how they describe a symbol, e.g. previous
# vs alias. Assume that HGNC annotations are correct, as HGNC is the
# naming consortium.  Remove entries where the value is the same
# (Regardless of how it is annotated)

ncbi_f <- ncbi_f %>%
    dplyr::anti_join(hgnc %>% dplyr::select(HGNC_ID, ENSEMBL_ID, value))

genes <- dplyr::bind_rows(hgnc, ncbi_f) %>%
    dplyr::arrange(HGNC_ID, ENSEMBL_ID, value) %>%
    dplyr::group_by(HGNC_ID) %>%
    tidyr::fill(UNIPROT_IDS, NCBI_ID, HGNC_NAME, .direction = "updown") %>%

    # Check if values are unambiguous
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_ID)) %>%
    dplyr::filter(n_genes == 1 | source == "HGNC") %>%
    dplyr::select(-n_genes)


# Check for disagreements between org_db and NCBI ----
ncbi_not_org_db <- ncbi_genes %>%
    dplyr::select(ENTREZ_ID, ENSEMBL_ID) %>%
    unique() %>%
    dplyr::anti_join(org_db %>%
                         dplyr::select(ENTREZ_ID, ENSEMBL_ID) %>%
                         unique())



# Genes in NCBI not org_db are where ENTREZ_ID is mapped to a novel gene or
# scaffold.



# Both org.db and ncbi can map one ENTREZ_ID to several ENSEMBL_IDs



# Add in the org.db genes -----




# Only want genes for which there is a HGNC / ENSEMBL combination in HGNC
org_f <- org_db %>%
    dplyr::semi_join(hgnc %>% dplyr::select(HGNC_ID, ENSEMBL_ID))

# NCBI and HGNC can differ in how they describe a symbol, e.g. previous
# vs alias. Assume that HGNC annotations are correct, as HGNC is the
# naming consortium.  Remove entries where the value is the same
# (Regardless of how it is annotated)
ncbi_f <- ncbi_f %>%
    dplyr::anti_join(hgnc %>% dplyr::select(HGNC_ID, ENSEMBL_ID, value))






















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



