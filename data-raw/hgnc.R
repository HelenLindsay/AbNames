# Download HGNC data ----
# Download gene groups and protein-coding genes tables from
# Human Genome Naming Consortium.  As at 01.06.2022, the groups
# table is not a subset of the protein-coding genes table.

# Example of a multi-subunit protein - LFA-1 from ITGAL and ITGB2

# HGNC proteins
hgnc_proteins_fname <- paste(c("http://ftp.ebi.ac.uk/pub/databases/",
                             "genenames/hgnc/tsv/locus_types/",
                             "gene_with_protein_product.txt"),
                           collapse = "")
hgnc_proteins_f <- tempfile()
download.file(hgnc_proteins_fname, destfile = hgnc_proteins_f)


# HGNC groups
hgnc_groups_fname <- paste(c("https://www.genenames.org/cgi-bin/genegroup/",
                             "download-all"), collapse = "")
hgnc_groups_f <- tempfile()
download.file(hgnc_groups_fname, destfile = hgnc_groups_f)

# Select relevant columns, rename and merge tables ----
hgnc_proteins <- readr::read_delim(hgnc_proteins_f)

hgnc_proteins <- hgnc_proteins[, c("hgnc_id", "symbol", "name", "alias_symbol",
                                   "prev_symbol", "ensembl_gene_id",
                                   "alias_name", "prev_name", "uniprot_ids")]

hgnc_proteins <- dplyr::rename(hgnc_proteins,
                               HGNC_ID = hgnc_id,
                               HGNC_NAME = name,
                               ENSEMBL_ID = ensembl_gene_id,
                               UNIPROT_ID = uniprot_ids,
                               HGNC_SYMBOL = symbol,
                               ALIAS = alias_symbol,
                               PREVIOUS_SYMBOL = prev_symbol,
                               ALIAS_NAME = alias_name,
                               PREVIOUS_NAME = prev_name)

hgnc_groups <- readr::read_delim(hgnc_groups_f)

hgnc_groups <- hgnc_groups %>%
    dplyr::rename(hgnc_groups,
                  HGNC_ID = `HGNC ID`,
                  HGNC_NAME = `Approved name`,
                  ENSEMBL_ID = `Ensembl gene ID`,
                  HGNC_SYMBOL = `Approved symbol`,
                  ALIAS = `Alias symbols`,
                  PREVIOUS_SYMBOL = `Previous symbols`,
                  ENTREZ_ID = `NCBI Gene ID`,
                  BIOTYPE = `Locus type`) %>%
    # Filter out pseudogenes and RNAs
    dplyr::filter(! grepl("^RNA|pseudogene|unknown", BIOTYPE)) %>%
    dplyr::select(hgnc_groups, -Status, -`BIOTYPE`, -`Chromosome`,
                             -`Vega gene ID`, `Group ID`, -`Group name`,
                             -`Group ID`)

hgnc_groups <-  %>%
    dplyr::mutate(across(c(PREVIOUS_SYMBOL, ALIAS), ~gsub(", ", "\\|", .x)))

hgnc <- dplyr::full_join(hgnc_proteins, hgnc_groups)


# Create and save long version of the HGNC table for querying ----

hgnc <- hgnc %>%
    dplyr::mutate(HGNC_SYMBOL2 = HGNC_SYMBOL) %>%
    tidyr::pivot_longer(c("HGNC_SYMBOL",
                          "HGNC_NAME",
                          "ALIAS",
                          "PREVIOUS_SYMBOL",
                          "ALIAS_NAME",
                          "PREVIOUS_NAME"),
                        names_to = "symbol_type") %>%
    dplyr::filter(! is.na(value)) %>%
    # Make one row per alias
    AbNames::splitUnnest(ab = "value", split = "\\|") %>%
    tidyr::unnest(cols = value) %>%
    dplyr::mutate(SOURCE = "HGNC") %>%
    dplyr::rename(HGNC_SYMBOL = HGNC_SYMBOL2) %>%
    unique() %>%

    # Remove entries that match the HGNC_SYMBOL, redundant
    dplyr::filter(! (value == HGNC_SYMBOL & ! symbol_type == "HGNC_SYMBOL")) %>%

    # Only keep ambiguous values if it's an official symbol
    dplyr::group_by(value) %>%
    dplyr::mutate(ngroups = n_distinct(HGNC_ID)) %>%
    dplyr::filter(ngroups == 1 | symbol_type == "HGNC_SYMBOL") %>%
    dplyr::select(-ngroups)


# TO DO:
# Remove genes that do not have a protein ID?


hgnc <- as.data.frame(hgnc)
usethis::use_data(hgnc, overwrite = TRUE, compress = "bzip2")
