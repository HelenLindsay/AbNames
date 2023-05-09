# To do: CD98, cellmarker only lists one gene, totalseq 2
# CELLMARKER CD98 should aggregate of SLC3A2 and SLC7A5

gene_aliases <- as_tibble(gene_aliases) %>%
    dplyr::filter(! SOURCE %in% c("CELLMARKER", "CSPA"))

# Cell marker database ---------------------------------------------

# Notes cellmarker v1

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

# Notes cellmarker v2
# Cellmarker v2.0 maps one marker name to one gene
# Entrez ID is sometimes incorrect, e.g. pointing to a pseudogene

# To do: make sure that cellmarker cellMarker col mapped to same Symbol/ID
# column is a correct, unambiguous alias
# Is ENTREZ symbol the HGNC symbol? - No, sometimes it's previous symbol
# or not found


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
                  SYMBOL = Symbol, # DON"T KNOW YET IF IT IS CORRECT
                  UNIPROT_ID = UNIPROTID) %>%
    dplyr::mutate(SOURCE = "CELLMARKER") %>%
    dplyr::mutate(UNIPROT_ID = na_if(UNIPROT_ID, "NA")) %>%
    # As CELLMARKER aggregates markers from other sources
    # trust primary sources first
    dplyr::filter(! value %in% gene_aliases$value)

# Use org.db to check gene type and remove entries with inconsistant
# entrez ids
org_db <- AnnotationDbi::select(hs,
                                keys = pull(cellmarker, ENTREZ_ID),
                                keytype = "ENTREZID",
                                columns = c("SYMBOL","GENETYPE")) %>%
    dplyr::rename(ENTREZ_ID = ENTREZID)



# CHECK IF ANTIGEN/SYMBOL/ENTREZ_ID COMB IS CONSISTENT






# Patch UNIPROT IDs ----
ga_patch <- gene_aliases %>%
    dplyr::select(any_of(colnames(cellmarker))) %>%
    dplyr::filter()



# Match by antigen for filling NA
ga_coalesce <- ga_patch[match(cellmarker$Antigen, ga_patch$Antigen), ]

# Remove Antigen for updating UNIPROT ID
ga_patch <- ga_patch %>%
    dplyr::select(-Antigen) %>%
    dplyr::filter(complete.cases(UNIPROT_ID, ENTREZ_ID)) %>%
    unique()

cellmarker <- cellmarker %>%
    # Fill NA values
    dplyr::coalesce(ga_coalesce) %>%

    # Filter again, some cellNames include an NA in a list of markers
    # (some aren't protein-coding, some would have an ENTREZ ID)
    dplyr::filter(if_all(c(ENTREZ_ID, UNIPROT_ID), ~! is.na(.x))) %>%
    dplyr::mutate(across(everything(), ~stringr::str_squish(.x))) %>%

    # Update UNIPROT IDs, matching by ENTREZ_ID
    dplyr::rows_update(ga_patch, by = "ENTREZ_ID", unmatched = "ignore") %>%

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
    dplyr::rename(HGNC_SYMBOL = ENTREZ_SYMBOL) %>%

    # Add HGNC_IDs and ENSEMBL_IDs
    dplyr::left_join(gene_aliases %>%
                         dplyr::select("HGNC_SYMBOL",
                                       "HGNC_ID",
                                       "ENSEMBL_ID") %>%
                         unique()) %>%

    # Re-aggregate the protein complexes
    unique() %>%
    dplyr::group_by(Antigen, cellName) %>%
    dplyr::mutate(across(c(HGNC_ID, HGNC_SYMBOL, ENSEMBL_ID, ENTREZ_ID,
                           proteinName, UNIPROT_ID),
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

# Add to aliases table, save
gene_aliases <- gene_aliases %>%
    dplyr::bind_rows(cellmarker)


# note still contains multiple antigen to same gene e.g. CD45 isoforms

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
#missing <- anti_join(cellmarker %>% dplyr::mutate(row_n = dplyr::row_number()),
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
    # No gene symbol
    dplyr::filter(! is.na(ENTREZ_SYMBOL)) %>%
    # Join by exact match to "value" in gene_aliases table
    # (either official symbol or alias)
    dplyr::inner_join(gene_aliases %>%
                          dplyr::select(HGNC_SYMBOL, HGNC_ID, ENSEMBL_ID,
                                        ENTREZ_ID, value) %>%
                          unique(),
                      by = c("ENTREZ_SYMBOL" = "value", "ENTREZ_ID")) %>%
    # Know from join that ENTREZ_SYMBOL must match an existing symbol or alias
    tidyr::pivot_longer(c(Antigen, UNIPROT_NAME), names_to = "symbol_type") %>%
    dplyr::filter(! (value %in% gene_aliases$value | is.na(value))) %>%
    dplyr::select(dplyr::any_of(colnames(gene_aliases)))


# Check for inconsistencies -----

# Some HGNC_ID / UNIPROT_ID combinations differ.
# Can be because gene_aliases UNIPROT_ID is NA, or because ID is out of date
aj <- cspa %>%
    dplyr::anti_join(gene_aliases, by = c("HGNC_ID", "ENTREZ_ID", "UNIPROT_ID"))


# Add to aliases table, save
gene_aliases <- gene_aliases %>%
    dplyr::bind_rows(cspa)

gene_aliases <- as.data.frame(gene_aliases)
usethis::use_data(gene_aliases, overwrite = TRUE, compress = "bzip2")
