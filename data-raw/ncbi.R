# Notes ----
#
# Entrez gene is a database for gene-specific information.  Includes curated
# and automatic integration of RefSeq (sequence-based annotation)
#
# NCBI sometimes maps one HGNC_ID to multiple Ensembl IDs.  In these cases,
# the chromosomal locations of the genes match.
# e.g. NAT-1 appears to be a multi-copy gene.
#
# Not all NCBI genes have a HGNC_SYMBOL or ENSEMBL_ID.
# Could these be processed pseudogenes?

# Generate long-format table from NCBI -----

library(tidyverse)
library(AbNames)
data(hgnc)

# Get NCBI genes and aliases

#fn <- file.path("https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens",
#                "Homo_sapiens.gene_info.gz")
#f <- tempfile()
#download.file(fn, destfile = f)

f <- "~/Analyses/CITEseq_curation/data/Homo_sapiens.gene_info.gz"
ncbi_genes <- readr::read_delim(f, na = c("", "NA", "-"))

# Keywords for filtering "Other" column
alias_grep <- "^CD|[Aa]ntigen|MHC|HLA|(T[- ]cell)|(B[- ]cell)|surface|immunoglo"


# Make long format NCBI table

ncbi_genes <- ncbi_genes %>%

    dplyr::filter(type_of_gene == "protein-coding",
                  if_all(.cols = c(Full_name_from_nomenclature_authority,
                                   Other_designations),
                         ~! grepl("pseudogene|uncharacterized|readthrough",
                                  .x))) %>%
    dplyr::rename(HGNC_SYMBOL = Symbol_from_nomenclature_authority,
                  HGNC_NAME = Full_name_from_nomenclature_authority,
                  NCBI_SYMBOL = Symbol,
                  ENTREZ_ID = GeneID,
                  ALIAS = Synonyms,
                  NCBI_NAME = description,
                  OTHER = Other_designations,
                  ENTREZ_ID = GeneID,
                  BIOTYPE = type_of_gene) %>%

    dplyr::mutate(HGNC_ID = stringr::str_extract_all(dbXrefs,
                                                    "(?<=HGNC:)HGNC:[0-9]+"),
                  ENSEMBL_ID = stringr::str_extract_all(dbXrefs,
                                                    "(?<=Ensembl:)ENSG[0-9]+"),
                  HGNC_ID = map_chr(HGNC_ID, AbNames:::.toString),
                  HGNC_ID = ifelse(HGNC_ID %in% c("NA", ""), NA, HGNC_ID),
                  ENTREZ_ID = as.character(ENTREZ_ID)) %>%

    # Assume that antibodies of interest are against proteins with HGNC IDs
    dplyr::filter(! is.na(HGNC_ID)) %>%

    tidyr::unnest(ENSEMBL_ID) %>%

    # Only keep genes for which there is a HGNC / ENSEMBL / ENTREZ comb in HGNC
    # (Note that not all HGNC IDs are in hgnc data set, e.g. non-protein-coding)

    # There are some differences in which ENSEMBL_ID is mapped to which HGNC_ID,
    # e.g. in HGNC HGNC:4883 -> ENSG00000000971 (official gene)
    #      in NCBI HGNC:4883 -> ENSG00000289697 (novel gene)
    dplyr::semi_join(hgnc %>% dplyr::select(HGNC_ID, HGNC_SYMBOL,
                                            ENSEMBL_ID, ENTREZ_ID)) %>%

    dplyr::mutate(# Only keep NCBI_SYMBOL if different from HGNC_SYMBOL
                  NCBI_SYMBOL = AbNames:::.noDups(NCBI_SYMBOL, HGNC_SYMBOL),
                  # Only keep NCBI_NAME if different from HGNC_NAME
                  NCBI_NAME = AbNames:::.noDups(NCBI_NAME, HGNC_NAME)) %>%

    # Only keep columns of interest
    dplyr::select(-`#tax_id`, -`map_location`, -`Modification_date`,
                  -Feature_type, -chromosome, -LocusTag,
                  -Nomenclature_status, -dbXrefs) %>%

    # Convert to long format, split aliases and other designations
    tidyr::pivot_longer(c(NCBI_SYMBOL, ALIAS, NCBI_NAME, OTHER),
                        names_to = "symbol_type") %>%
    dplyr::filter(! is.na(value)) %>%
    AbNames::splitUnnest(ab = "value", "\\|") %>%

    # Filter "Other" column for terms related to surface markers
    dplyr::filter(symbol_type == "OTHER" & grepl(alias_grep, value) |
                      ! symbol_type == "OTHER",
                  # Remove entries such as "antigen-like", "immunoglobulin like"
                  ! (symbol_type == "OTHER" & grepl("like", value))) %>%

    dplyr::mutate(SOURCE = "NCBI") %>%

    # Remove aliases that map to more than one HGNC_SYMBOL
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_SYMBOL)) %>%
    dplyr::filter(n_genes == 1) %>%
    dplyr::select(-n_genes)


ncbi_genes <- as.data.frame(ncbi_genes)
usethis::use_data(ncbi_genes, overwrite = TRUE, compress = "bzip2")
