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

# Setup -----

library(tidyverse)
library(readxl)

source("ncbi.R")
source("org_db.R")

# --------------------------------------------------------------------
# Make HGNC/ENSEBML IDs from HGNC consistent with Ensembl -----

# Two cases:
# - ENSEMBL / HGNC disagree -> use ENSEMBL ID from ENSEMBL
# - ENSEMBL has additional matches -> keep all protein-coding

# Manually checked 15/06/22, usually disagreements are because
# HGNC ENSEMBL_ID is obsolete

bm_ids <- dplyr::select(bm, ENSEMBL_ID, HGNC_ID, HGNC_SYMBOL, BIOTYPE) %>%
    unique()

# 28 genes differ in HGNC_ID / ENSEMBL_ID combination
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

# These are the genes where HGNC/ENSEMBL disagree
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

# Add both copies into hgnc
temp <- hgnc %>% dplyr::semi_join(multi_gene)
multi_gene_patch <- multi_gene %>% dplyr::anti_join(temp)
temp <- dplyr::rows_update(temp,
                           multi_gene_patch %>% dplyr::select(-BIOTYPE),
                           by = c("HGNC_ID", "HGNC_SYMBOL"))
hgnc <- hgnc %>% bind_rows(temp)

# Note that the mappings to non-protein-coding genes haven't been removed
# e.g. ENSG00000271672 transcribed_processed_pseudogene

# -----------------------------------------------------------------------
# Add Entrez IDs -----

prep_merge <- function(df, hgnc){
    df %>%
        dplyr::select(ENSEMBL_ID, ENTREZ_ID, HGNC_SYMBOL) %>%
        dplyr::semi_join(hgnc) %>%
        # Want Ensembl, Entrez and Symbol for merging
        na.omit() %>%
        unique()
}

common_or_missing <- function(x, y){
    sj <- dplyr::semi_join(x, y)
    aj <- dplyr::anti_join(x, y)
    # Are the entries in the anti-join inconsistent or missing?

    # Find entries in y where any column value is in the anti_join
    uj <- union_join(y, aj)

    # Find entries in the anti-join where any entry is in y (using above)
    to_rm <- union_join(aj, uj)

    aj <- anti_join(aj, to_rm)
    return(dplyr::bind_rows(sj, aj))
}


bm_x <- prep_merge(bm, hgnc)
ncbi_x <- prep_merge(ncbi_genes, hgnc)
org_db_x <- prep_merge(org_db, hgnc)

bm_ncbi <- common_or_missing(bm_x, ncbi_x)
bm_orgdb <- common_or_missing(bm_x, org_db_x)
ncbi_bm <- common_or_missing(ncbi_x, bm_x)
ncbi_orgdb <- common_or_missing(ncbi_x, org_db_x)
orgdb_bm <- common_or_missing(org_db_x, bm_x)
orgdb_ncbi <- common_or_missing(org_db_x, ncbi_x)

consistent <- Reduce(dplyr::full_join,
                     list(bm_ncbi, bm_orgdb, ncbi_bm, ncbi_orgdb,
                                  orgdb_bm, orgdb_ncbi))

hgnc <- hgnc %>% dplyr::left_join(consistent)


# ---------------------------------------------------------------------------
# Select the aliases that aren't already present in hgnc -----

bm_novel <- bm %>%
    dplyr::mutate(symbol_type = "ALIAS") %>%
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

# --------------------------------------------------------------------------

# Add novel aliases -----

# Add novel aliases from biomaRt and NCBI
# (Note - have previously checked that HGNC_ID / ENSEMBL_ID combinations
# are correct)
hgnc <- hgnc %>%
    dplyr::bind_rows(bm_novel, ncbi_novel, org_db_novel) %>%

    # Fill Entrez IDs from new aliases
    dplyr::group_by(HGNC_ID) %>%
    tidyr::fill(ENTREZ_ID, HGNC_NAME, UNIPROT_IDS, .direction = "updown") %>%

    # Fill HGNC IDs (missing from org_db)
    dplyr::group_by(ENSEMBL_ID) %>%
    tidyr::fill(HGNC_ID, HGNC_NAME, UNIPROT_IDS, .direction = "updown") %>%

    # New values may come from more than one source, aggregate
    dplyr::group_by(HGNC_ID, ENSEMBL_ID, UNIPROT_IDS, HGNC_SYMBOL, ENTREZ_ID,
                    symbol_type, value) %>%
    dplyr::summarise(SOURCE = toString(SOURCE), .groups = "keep") %>%

    # Check if values are unambiguous
    # (newly added aliases may create new ambiguities, only 3 are ambiguous)
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_ID)) %>%
    dplyr::filter(n_genes == 1 | SOURCE == "HGNC") %>%
    dplyr::select(-n_genes)


#> table(is.na(hgnc$ENTREZ_ID)) # Note - not the number of genes because aliases
# 122569  36302

