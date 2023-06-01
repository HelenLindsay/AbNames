# Download HGNC data ----
# Download gene groups and protein-coding genes tables from
# Human Genome Naming Consortium.  As at 01.06.2022, the groups
# table is not a subset of the protein-coding genes table.

# Example of a multi-subunit protein - LFA-1 from ITGAL and ITGB2
existing <- ls()

# HGNC proteins
hgnc_proteins_fname <- paste0("http://ftp.ebi.ac.uk/pub/databases/",
                             "genenames/hgnc/tsv/locus_types/",
                             "gene_with_protein_product.txt")
hgnc_proteins_f <- sprintf("%s/hgnc_gene_with_protein_product_%s.txt",
                           downloads, Sys.Date())
if (! file.exists(hgnc_proteins_f)){
    download.file(hgnc_proteins_fname, destfile = hgnc_proteins_f)
}

# HGNC groups
hgnc_groups_fname <- paste(c("https://www.genenames.org/cgi-bin/genegroup/",
                             "download-all"), collapse = "")
hgnc_groups_f <- sprintf("%s/hgnc_all_groups_%s.csv", downloads, Sys.Date())
if (! file.exists(hgnc_groups_f)){
    download.file(hgnc_groups_fname, destfile = hgnc_groups_f)
}

# Select relevant columns, rename and merge tables ----
hgnc_proteins <- readr::read_delim(hgnc_proteins_f)

hgnc_proteins <- hgnc_proteins[, c("hgnc_id", "symbol", "name", "alias_symbol",
                                   "prev_symbol", "ensembl_gene_id",
                                   "entrez_id", "alias_name", "prev_name",
                                   "uniprot_ids", "locus_type")]

hgnc_proteins <- dplyr::rename(hgnc_proteins,
                               HGNC_ID = hgnc_id,
                               HGNC_NAME = name,
                               ENSEMBL_ID = ensembl_gene_id,
                               ENTREZ_ID = entrez_id,
                               UNIPROT_ID = uniprot_ids,
                               HGNC_SYMBOL = symbol,
                               ALIAS = alias_symbol,
                               PREVIOUS_SYMBOL = prev_symbol,
                               ALIAS_NAME = alias_name,
                               PREVIOUS_NAME = prev_name,
                               BIOTYPE = locus_type)

hgnc_groups <- readr::read_delim(hgnc_groups_f)

hgnc_groups <- hgnc_groups %>%
    dplyr::rename(HGNC_ID = `HGNC ID`,
                  HGNC_NAME = `Approved name`,
                  ENSEMBL_ID = `Ensembl gene ID`,
                  HGNC_SYMBOL = `Approved symbol`,
                  ALIAS = `Alias symbols`,
                  PREVIOUS_SYMBOL = `Previous symbols`,
                  ENTREZ_ID = `NCBI Gene ID`,
                  BIOTYPE = `Locus type`) %>%
    # Filter out pseudogenes and RNAs
    dplyr::filter(! grepl("^RNA|pseudogene|unknown|retrovirus|readthrough",
                          BIOTYPE),
                  ! grepl("^MT-", HGNC_SYMBOL)) %>%
    dplyr::select(-Status, -`Chromosome`, -`Vega gene ID`,
                  `Group ID`, -`Group name`, -`Group ID`) %>%
    dplyr::mutate(across(c(PREVIOUS_SYMBOL, ALIAS), ~gsub(", ", "\\|", .x)))

hgnc <- dplyr::full_join(hgnc_proteins, hgnc_groups)

# Check that one HGNC ID maps to one ENTREZ/ENSEMBL ID ----
temp <- hgnc %>%
    AbNames:::nPerGroup(group = "HGNC_ID",
                       col = c("ENTREZ_ID", "ENSEMBL_ID")) %>%
    dplyr::filter(nENTREZ_ID > 1 | nENSEMBL_ID > 1)

stopifnot(nrow(temp) == 0)


# Check that one HGNC symbol maps to one ID ----
temp <- hgnc %>%
    AbNames:::nPerGroup(group = "HGNC_SYMBOL",
                        col = c("HGNC_ID", "ENTREZ_ID", "ENSEMBL_ID")) %>%
    dplyr::filter(nENTREZ_ID > 1 | nENSEMBL_ID > 1 | nHGNC_ID > 1)

stopifnot(nrow(temp) == 0)


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
    dplyr::mutate(SOURCE = "HGNC",
                  BIOTYPE = ifelse(BIOTYPE == "gene with protein product",
                                   "protein_coding", BIOTYPE)) %>%
    dplyr::rename(HGNC_SYMBOL = HGNC_SYMBOL2) %>%
    unique() %>%

    # Remove entries that match the HGNC_SYMBOL, redundant
    dplyr::filter(! (value == HGNC_SYMBOL & ! symbol_type == "HGNC_SYMBOL")) %>%

    # Only keep ambiguous values if it's an official symbol
    dplyr::group_by(value) %>%
    dplyr::mutate(ngroups = n_distinct(HGNC_ID)) %>%
    dplyr::filter(ngroups == 1 | symbol_type == "HGNC_SYMBOL") %>%
    dplyr::select(-ngroups) %>%

    # Make ENTREZ ID a character for joining with other tables
    dplyr::mutate(ENTREZ_ID = as.character(ENTREZ_ID)) %>%
    dplyr::ungroup()

write_csv(hgnc, file = sprintf("%s/hgnc.csv", downloads))
rm(list = setdiff(ls(), c(existing, "hgnc", "existing")))

#hgnc <- as.data.frame(hgnc)
#usethis::use_data(hgnc, overwrite = TRUE, compress = "bzip2")

# Possible changes:
# Remove genes that do not have a protein ID?
# Filter ambiguous previous symbols then recheck for ambiguity?
