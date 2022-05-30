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
    dplyr::filter(type_of_gene == "protein-coding") %>%
      dplyr::rename(HGNC_SYMBOL = Symbol_from_nomenclature_authority,
                    HGNC_NAME = Full_name_from_nomenclature_authority,
                    NCBI_SYMBOL = Symbol,
                    ALIAS = Synonyms,
                    NCBI_DESC = description,
                    EXT_DESC = Other_designations) %>%
    dplyr::mutate(HGNC_ID = stringr::str_extract_all(dbXrefs,
                                                    "(?<=HGNC:)HGNC:[0-9]+"),
                  ENSEMBL_ID = stringr::str_extract_all(dbXrefs,
                                                    "(?<=Ensembl:)ENSG[0-9]+"),
                  across(c(HGNC_ID, ENSEMBL_ID),  ~ map_chr(.x, toString)),
                  NCBI_SYMBOL = ifelse(NCBI_SYMBOL == HGNC_SYMBOL, NA,
                                       NCBI_SYMBOL)) %>%
    dplyr::select(-`#tax_id`, -`map_location`, -`type_of_gene`,
                  -`Modification_date`, -Feature_type, -chromosome,
                  -LocusTag, -Nomenclature_status, -dbXrefs, -GeneID) %>%

    # Remove the descriptions?
    tidyr::pivot_longer(c(HGNC_SYMBOL, NCBI_SYMBOL, ALIAS,
                          NCBI_DESC, EXT_DESC), names_to = "symbol_type") %>%
    dplyr::filter(! is.na(value)) %>%
    dplyr::mutate(value = strsplit(value, "\\|")) %>%
    tidyr::unnest(value)



