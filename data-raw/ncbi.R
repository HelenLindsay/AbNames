# Notes ----
#
# NCBI sometimes maps one HGNC_ID to multiple Ensembl IDs.  In these cases,
# the chromosomal locations of the genes match.
# e.g. NAT-1 appears to be a multi-copy gene.
#
#library(org.Hs.eg.db)
#org.db <- org.Hs.eg.db
#ensembl_locs <- AnnotationDbi::select(org.db,
#                                      keys = ncbi_genes$ENSEMBL_ID,
#                                      keytype = "ENSEMBL",
#                                      columns = c("MAP", "SYMBOL")) %>%
#    dplyr::rename(ENSEMBL_ID = ENSEMBL,
#                  ENTREZ_SYMBOL = SYMBOL)
#
# Not all NCBI genes have a HGNC_SYMBOL or ENSEMBL_ID.
# Could these be processed pseudogenes?

# Generate long-format table from NCBI -----

library(tidyverse)
library(AbNames)

# Get NCBI genes and aliases

fn <- file.path("https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens",
                "Homo_sapiens.gene_info.gz")
f <- tempfile()
download.file(fn, destfile = f)

ncbi_genes <- readr::read_delim(f, na = c("", "NA", "-"))

# Make long format NCBI table

ncbi_genes <- ncbi_genes %>%
    dplyr::filter(type_of_gene == "protein-coding",
                  ! grepl("pseudogene|uncharacterized|readthrough",
                          Other_designations)) %>%
      dplyr::rename(HGNC_SYMBOL = Symbol_from_nomenclature_authority,
                    HGNC_NAME = Full_name_from_nomenclature_authority,
                    NCBI_SYMBOL = Symbol,
                    NCBI_ID = GeneID,
                    ALIAS = Synonyms,
                    NCBI_NAME = description,
                    OTHER = Other_designations) %>%
    dplyr::mutate(HGNC_ID = stringr::str_extract_all(dbXrefs,
                                                    "(?<=HGNC:)HGNC:[0-9]+"),
                  ENSEMBL_ID = stringr::str_extract_all(dbXrefs,
                                                    "(?<=Ensembl:)ENSG[0-9]+"),
                  HGNC_ID = map_chr(HGNC_ID, toString),
                  HGNC_ID = ifelse(HGNC_ID %in% c("NA", ""), NA, HGNC_ID)) %>%
    # Assume that antibodies of interest are against proteins with HGNC IDs
    dplyr::filter(! is.na(HGNC_ID)) %>%
    tidyr::unnest(ENSEMBL_ID) %>%
    dplyr::mutate(# Only keep NCBI_SYMBOL if different from HGNC_SYMBOL
                  NCBI_SYMBOL = AbNames:::.noDups(NCBI_SYMBOL, HGNC_SYMBOL),
                  # Only keep NCBI_NAME if different from HGNC_NAME
                  NCBI_NAME = AbNames:::.noDups(NCBI_NAME, HGNC_NAME)) %>%
    dplyr::select(-`#tax_id`, -`map_location`, -`type_of_gene`,
                  -`Modification_date`, -Feature_type, -chromosome, -LocusTag,
                  -Nomenclature_status, -dbXrefs) %>%
    tidyr::pivot_longer(c(NCBI_SYMBOL, ALIAS, NCBI_NAME, OTHER),
                        names_to = "symbol_type") %>%
    dplyr::filter(! is.na(value)) %>%
    AbNames::splitUnnest(ab = "value", "\\|") %>%
    # Only keep "Other" column if it refers to an antigen or a CD molecule
    dplyr::filter(symbol_type == "OTHER" & grepl("[Aa]ntigen|^CD", value) |
                      ! symbol_type == "OTHER") %>%
    dplyr::mutate(source = "NCBI") %>%
    # Remove aliases that map to more than one HGNC_SYMBOL
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_SYMBOL)) %>%
    dplyr::filter(n_genes == 1) %>%
    dplyr::select(-n_genes)
