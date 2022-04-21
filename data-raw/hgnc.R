# Download HGNC data ----
# Download gene groups and protein-coding genes tables from
# Human Genome Naming Consortium.  As at 21.04.2022, the groups
# table is not a subset of the protein-coding genes table.

# HGNC proteins
hgnc_proteins_fname <- paste(c("http://ftp.ebi.ac.uk/pub/databases/",
                             "genenames/hgnc/tsv/locus_types/",
                             "gene_with_protein_product.txt"),
                           collapse = "")
hgnc_proteins <- tempfile()
download.file(hgnc_proteins_fname, destfile = hgnc_proteins)


# HGNC groups
hgnc_groups_fname <- paste(c("https://www.genenames.org/cgi-bin/genegroup/",
                             "download-all"), collapse = "")
hgnc_groups <- tempfile()
download.file(hgnc_groups_fname, destfile = hgnc_groups)

# Select relevant columns, rename and merge tables ----
hgnc_proteins <- readr::read_delim(hgnc_proteins)

hgnc_proteins <- hgnc_proteins[, c("hgnc_id", "symbol", "name", "alias_symbol",
                                   "prev_symbol", "ensembl_gene_id", "location",
                                   "alias_name", "prev_name", "uniprot_ids")]

hgnc_proteins <- dplyr::rename(hgnc_proteins,
                               HGNC_ID = hgnc_id,
                               HGNC_NAME = name,
                               ENSEMBL_ID = ensembl_gene_id,
                               UNIPROT_IDS = uniprot_ids,
                               HGNC_SYMBOL = symbol,
                               Chromosome = location)
hgnc_proteins <- dplyr::mutate(hgnc_proteins,
                               alias_symbol = gsub("\\|", ", ", alias_symbol),
                               prev_symbol = gsub("\\|", ", ", prev_symbol))

hgnc_groups <- readr::read_delim(hgnc_groups)

hgnc_groups <- dplyr::select(hgnc_groups, -Status, -`Locus type`,
                             -`NCBI Gene ID`, -`Vega gene ID`,
                             `Group ID`, -`Group name`, -`Group ID`)

hgnc_groups <- dplyr::rename(hgnc_groups,
                             HGNC_ID = `HGNC ID`,
                             HGNC_NAME = `Approved name`,
                             ENSEMBL_ID = `Ensembl gene ID`,
                             HGNC_SYMBOL = `Approved symbol`,
                             alias_symbol = `Alias symbols`,
                             prev_symbol = `Previous symbols`)

hgnc <- dplyr::full_join(hgnc_proteins, hgnc_groups)
hgnc <- unique(hgnc)
hgnc <- data.frame(hgnc)
usethis::use_data(hgnc, overwrite = TRUE)
