# Notes -----

# Genes in NCBI not org_db are where ENTREZ_ID is mapped to a novel gene or
# scaffold.

# Both org.db and ncbi can map one ENTREZ_ID to several ENSEMBL_IDs

# NCBI contains Entrez IDs not in org_db and vice versa,
# usually when they differ it is for this reason

# Example of gene with ambiguous aliases -
# CCL4L1/2 both have alias AT744.2

# NCBI genes has more cases where one HGNC symbol is mapped to several
# ENSEMBL IDs

# Not all HGNC Ensembl IDs are in the biomaRt table
# (some are not protein-coding)

# ENSG00000203668: Entrez gene symbol CHML, HGNC symbol OPN3 -
# both are valid HGNC symbols

# Within hgnc, one Uniprot ID can map to multiple HGNC IDs, e.g. Q9GZY0

# TO DO: FILTER OR FILL MISSING HGNC ENSEMBL_IDS

# Setup -----

library(tidyverse)
library(readxl)

source("ncbi.R")
source("entrez_ensembl.Rs")

# TO DO: what happened to biotype?
# --------------------------------------------------------------------
# Make HGNC/ENSEMBL IDs from HGNC consistent with Ensembl -----

# Two cases:
# - ENSEMBL / HGNC disagree -> use ENSEMBL ID from ENSEMBL
# - ENSEMBL has additional matches -> keep all protein-coding

# Manually checked 15/06/22, usually disagreements are because
# HGNC ENSEMBL_ID is obsolete

bm_ids <- dplyr::select(bm, ENSEMBL_ID, HGNC_ID, HGNC_SYMBOL, BIOTYPE) %>%
    unique()

# 15 genes differ in HGNC_ID / ENSEMBL_ID combination
bm_ens_diff <- dplyr::anti_join(bm_ids, hgnc_ids,
                                by = c("HGNC_ID", "ENSEMBL_ID"))

# Get the Ensembl IDs that HGNC assigns
hgnc_ens_diff <- hgnc_ids %>%
    union_join(bm_ens_diff %>% dplyr::select(-BIOTYPE))

# Check if these are also in the biomart table
# (2 ENSEMBL_IDs to 1 HGNC_SYMBOL)
multi_gene <- bm_ids %>%
    # HGNC_ID/ENSEMBL_ID agrees with HGNC
    dplyr::semi_join(hgnc_ens_diff)

# These are protein-coding genes where HGNC/ENSEMBL disagree (not multi-gene)
ens_patch <- bm_ens_diff %>%
    dplyr::filter(! HGNC_ID %in% multi_gene$HGNC_ID) %>%
    dplyr::filter(BIOTYPE == "protein_coding") %>%
    dplyr::select(-BIOTYPE)

# Overwrite the Ensembl IDs in HGNC if HGNC/ENSEMBL disagree
hgnc <- hgnc %>%
    dplyr::rows_update(ens_patch, by = c("HGNC_ID", "HGNC_SYMBOL"))

# Get both members where Ensembl gives two mappings to the same HGNC symbol
multi_gene <- bm_ids %>%
    dplyr::filter(HGNC_ID %in% multi_gene$HGNC_ID)

# Create a data frame with the alternative ENSEMBL IDs and add both
# copies into hgnc
temp <- hgnc %>%
    dplyr::semi_join(multi_gene,
                     by = c("ENSEMBL_ID", "HGNC_ID", "HGNC_SYMBOL")) %>%
    dplyr::select(-BIOTYPE, -ENSEMBL_ID)

hgnc <- hgnc %>% anti_join(temp)

temp <- full_join(temp, multi_gene)

hgnc <- hgnc %>% bind_rows(temp)


# ---------------------------------------------------------------------------
# Select the aliases that aren't already present in hgnc -----

bm_novel <- bm %>%
    dplyr::rename(value = ALIAS) %>%
    dplyr::filter(! value == "") %>%
    dplyr::anti_join(hgnc, by = c("HGNC_SYMBOL", "value")) %>%
    # Only keep if HGNC_SYMBOL/ENSEMBL_ID combinations are correct
    dplyr::semi_join(hgnc, by = c("HGNC_SYMBOL", "ENSEMBL_ID"))

org_db_novel <- org_db %>%
    dplyr::mutate(symbol_type = "ALIAS") %>%
    dplyr::anti_join(hgnc, by = c("HGNC_SYMBOL", "value")) %>%
    # Only want genes for which there is a HGNC / ENSEMBL combination in HGNC
    dplyr::semi_join(hgnc %>% dplyr::select(HGNC_ID, ENSEMBL_ID))

# NCBI and HGNC can differ in how they describe a symbol, e.g. previous
# vs alias. Assume that HGNC annotations are correct, as HGNC is the
# naming consortium.  Remove entries where the value is the same
# (Regardless of what type of symbol it is considered to be)
ncbi_novel <- ncbi_genes %>%
    dplyr::anti_join(hgnc %>% dplyr::select(HGNC_ID, ENSEMBL_ID, value))

# --------------------------------------------------------------------------

# Add novel aliases -----

# Add novel aliases from biomaRt and NCBI
# (Note - have previously checked that HGNC_ID / ENSEMBL_ID combinations
# are correct)
hgnc <- hgnc %>%
    dplyr::bind_rows(bm_novel, ncbi_novel, org_db_novel) %>%

    # Fill Entrez IDs from new aliases
    dplyr::group_by(HGNC_ID) %>%
    tidyr::fill(ENTREZ_ID, HGNC_NAME, UNIPROT_ID, .direction = "updown") %>%

    # Fill HGNC IDs (missing from org_db)
    dplyr::group_by(ENSEMBL_ID) %>%
    tidyr::fill(HGNC_ID, HGNC_NAME, UNIPROT_ID, .direction = "updown") %>%

    # New values may come from more than one source, aggregate
    dplyr::group_by(HGNC_ID, ENSEMBL_ID, UNIPROT_ID, HGNC_SYMBOL, ENTREZ_ID,
                    symbol_type, value) %>%
    dplyr::summarise(SOURCE = toString(SOURCE), .groups = "keep") %>%

    # Check if values are unambiguous
    # (newly added aliases may create new ambiguities, only 3 are ambiguous)
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_ID)) %>%
    dplyr::filter(n_genes == 1 | SOURCE == "HGNC") %>%
    dplyr::select(-n_genes) %>%

    # Aggregate the protein IDs
    dplyr::group_by(HGNC_ID, HGNC_SYMBOL, ENSEMBL_ID, ENTREZ_ID) %>%
    dplyr::mutate(UNIPROT_ID = paste0(sort(unique(UNIPROT_ID)),
                                         collapse = "|")) %>%
    ungroup() %>%
    unique()

gene_aliases <- as.data.frame(hgnc)
usethis::use_data(gene_aliases, overwrite = TRUE, compress = "bzip2")
