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


# Cell marker database ---------------------------------------------

# http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.jsp

gsubCellmarker <- function(x){
    p1_matches <- stringr::str_extract_all(x, "\\[([^\\[]+)\\]")
    p1_sub <- gsub(", ", "\\|", unlist(p1_matches))
    p1_sub <- gsub("\\[|\\]", "", p1_sub)
    p1_sub <- relist(p1_sub, p1_matches)

    x <- cellmarker$geneSymbol
    m <- gregexpr("\\[([^\\[]+)\\]", x[60:65])
    regmatches(x[60:65], m) <- p1_sub[60:65]
    return(x)
}

# Cellmarker
cellmarker_loc <- "http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt"
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
                           proteinName, proteinID))

protein_complexes <- cellmarker %>%
    dplyr::filter(grepl("\\|", geneSymbol)) %>%
    dplyr::select(cellMarker, geneSymbol, geneID, proteinName, proteinID) %>%
    unique()


# Row Epithelial cell starting Adhesion molecules -
# Number of cell markers is not equal to number of gene (groups),
# is there an extra , Adhesion molecules, LFA1, Adhesion molecules LFA2?


#dplyr::rename(ENTREZ_ID = geneID,
#             UNIPROT_ID = proteinID) %>%


# Cell surface protein atlas ---------------------------------------

# http://wlab.ethz.ch/cspa/#abstract

library(readxl)

# Table of validated surfaceome proteins
cspa_loc <- "http://wlab.ethz.ch/cspa/data/S2_File.xlsx"
cspa_fname <- "~/Analyses/CITEseq_curation/data/cspa.xlsx"
download.file(cspa_loc, destfile = cspa_fname)

# Sheet A has human proteins and entrez IDs


# EBI complex portal ----------------------------------------------

# Human complex tab
ebi <- paste("http://ftp.ebi.ac.uk/pub/databases/intact/complex/",
             "current/complextab/9606.tsv")




