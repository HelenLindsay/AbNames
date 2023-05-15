# CD16 alias for CD16a

gene_aliases <- as_tibble(gene_aliases) %>%
    dplyr::filter(! SOURCE %in% c("CELLMARKER", "CSPA"))

# Cell marker database ---------------------------------------------

# Notes cellmarker v1 ----

# Note: Isoform name may be mapped to gene name, e.g. CD45RO -> PTPRC
# Some antigens have no annotation information, e.g. CD45.1
# Table contains annotation errors, e.g. row shift error in UNIPROT_IDs

# Row Epithelial cell starting Adhesion molecules -
# Number of cell markers is not equal to number of gene (groups),
# is there an extra , Adhesion molecules, LFA1, Adhesion molecules LFA2?

# V1 There are some errors in the uniprot IDs, e.g. P08920 for CD2 is a mouse ID

# Cellmarker 1.0 is incorrect for CD77 / A4GALT
# A4GALT synthesises CD77, CD77 is not a protein
#cellmarker_exclude <- c("CD77")

# Cellmarker v1
#cellmarker_loc <- paste0("http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/",
#                         "Human_cell_markers.txt")

# Notes cellmarker v2 ----
# Symbol / Entrez wrong from CD16
# value / Symbol wrong for SEP9 / spt8

# Cellmarker v2.0 maps one marker name to one gene
# Entrez ID is sometimes incorrect, e.g. pointing to a pseudogene

# To do: make sure that cellmarker cellMarker col mapped to same Symbol/ID
# column is a correct, unambiguous alias
# Is ENTREZ symbol the HGNC symbol? - No, sometimes it's previous symbol
# or not found
# To do: CD98, cellmarker only lists one gene, totalseq 2
# CELLMARKER CD98 should aggregate of SLC3A2 and SLC7A5

# Cellmarker 2.0
cellmarker_loc <- paste0("http://bio-bigdata.hrbmu.edu.cn/CellMarker/",
                         "CellMarker_download_files/file/Cell_marker_Human.xlsx")
cellmarker_fname <- sprintf("inst/extdata/CellMarker_human_%s.xlsx", Sys.Date())
download.file(cellmarker_loc, destfile = cellmarker_fname)

cellmarker <- readxl::read_xlsx(cellmarker_fname) %>%
    dplyr::filter(Genetype == "protein_coding") %>%
    dplyr::select(marker, Symbol, GeneID, UNIPROTID) %>%
    unique() %>%
    dplyr::rename(value = marker,
                  ENTREZ_ID = GeneID,
                  HGNC_SYMBOL = Symbol, # Don't know if it is correct yet
                  UNIPROT_ID = UNIPROTID) %>%
    dplyr::mutate(SOURCE = "CELLMARKER",
                  symbol_type = "ALIAS",
                  BIOTYPE = "protein-coding") %>%
    dplyr::mutate(UNIPROT_ID = na_if(UNIPROT_ID, "NA")) %>%
    # As CELLMARKER aggregates markers from other sources
    # trust primary sources first
    dplyr::filter(! value %in% gene_aliases$value) %>%
    dplyr::semi_join(gene_aliases, by = c("HGNC_SYMBOL", "ENTREZ_ID"))

# Add HGNC_ID and ENSEML_ID
ga_patch <- gene_aliases %>%
    dplyr::select(HGNC_ID, ENSEMBL_ID, HGNC_SYMBOL, ENTREZ_ID) %>%
    unique()

cellmarker <- cellmarker %>%
    dplyr::left_join(ga_patch, relationship = "many-to-one") %>%
    dplyr::filter(! is.na(HGNC_ID))

# Add to aliases table, save
gene_aliases <- gene_aliases %>%
    dplyr::bind_rows(cellmarker)

# Cell surface protein atlas ---------------------------------------

# http://wlab.ethz.ch/cspa/#abstract

# Notes:
# Putative Ig-like domain-containing protein / 374383
# ODZ3 / 55714 = TENM3 checked manually via gene cards

# Most entries can be matched to HGNC via UNIPROT_ID / HGNC_SYMBOL
# Those that can't include HLA and T-cell receptors

# Does not appear to contain protein complex names

# Table of validated surfaceome proteins

cspa_loc <- "http://wlab.ethz.ch/cspa/data/S2_File.xlsx"
cspa_fname <- sprintf("inst/extdata/cspa_%s.xlsx", Sys.Date())
download.file(cspa_loc, destfile = cspa_fname)

# Sheet A has human proteins and entrez IDs
cspa <- readxl::read_xlsx(cspa_fname) %>%
    dplyr::select(UP_Protein_name, CD, ENTREZ_gene_ID,
                  `ENTREZ gene symbol`, ID_link) %>%
    dplyr::rename(UNIPROT_NAME = UP_Protein_name,
                  UNIPROT_ID = ID_link,
                  ENTREZ_ID = ENTREZ_gene_ID,
                  ENTREZ_SYMBOL = `ENTREZ gene symbol`,
                  Antigen = CD) %>%
    dplyr::mutate(SOURCE = "CSPA",
                  Antigen = na_if(Antigen, "no"),
                  ENTREZ_ID = as.character(ENTREZ_ID),
                  across(c(ENTREZ_SYMBOL, ENTREZ_ID), ~na_if(.x, "0"))) %>%
    # No gene symbol
    dplyr::filter(! is.na(ENTREZ_SYMBOL)) %>%

    # Join by exact match to "value" in gene_aliases table
    # (either official symbol or alias)
    dplyr::inner_join(gene_aliases %>%
                          dplyr::select(HGNC_SYMBOL, HGNC_ID, ENSEMBL_ID,
                                        ENTREZ_ID, value) %>%
                          unique(),
                      by = c("ENTREZ_SYMBOL" = "value", "ENTREZ_ID"),
                      relationship = "many-to-many") %>%

    # If the ENTREZ_SYMBOL is not the official symbol, list as an alias
    dplyr::mutate(ENTREZ_SYMBOL =
                      ifelse(ENTREZ_SYMBOL == HGNC_SYMBOL,
                             NA, ENTREZ_SYMBOL)) %>%

    # Know from join that ENTREZ_SYMBOL must match an existing symbol or alias
    tidyr::pivot_longer(c(Antigen, UNIPROT_NAME, ENTREZ_SYMBOL),
                        names_to = "symbol_type") %>%
    dplyr::filter(! (value %in% gene_aliases$value | is.na(value)))

# Check for inconsistencies -----

# Some HGNC_ID / UNIPROT_ID combinations differ.
# Can be because gene_aliases UNIPROT_ID is NA, or because ID is out of date
# (all agree on HGNC_SYMBOL and ENTREZ_ID by construction)
aj <- cspa %>%
    dplyr::anti_join(gene_aliases, by = c("HGNC_ID", "ENTREZ_ID", "UNIPROT_ID"))

# Add to aliases table, save
gene_aliases <- gene_aliases %>%
    dplyr::bind_rows(cspa)

write_csv(as.data.frame(hgnc), file = "inst/extdata/gene_aliases.csv")
rm(list = setdiff(ls(), c("existing","gene_aliases")))

#gene_aliases <- as.data.frame(gene_aliases)
#usethis::use_data(gene_aliases, overwrite = TRUE, compress = "bzip2")
