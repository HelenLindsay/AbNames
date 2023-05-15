# CD97 matches ADRE5 / HGNC:1711

# Notes -----

# PROBLEMS:
# THERE ARE DUPLICATED VALUES WITH SOURCE NCBI, ORGDB / ORGDB, e.g. KRAS
# HGNC:17134 - HGNC has a UNIPROT ID, why is there no source == HGNC?
# ARMC9 / ENSG00000135931 for isotype control
# Antigen with e.g. CD3.1 - not safe to remove dot if both sides are number

# Ensembl can map same gene to multiple HGNC symbols, e.g. ENSG00000276085
# to CCL3L1 / CCL3L3.  HGNC maps to different ENSEMBL IDs,

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

existing = ls()

hgnc <- dplyr::mutate(hgnc, ENTREZ_ID = as.character(ENTREZ_ID))
# --------------------------------------------------------------------
# Make HGNC/ENSEMBL IDs from HGNC consistent with Ensembl -----

# Two cases:
# - ENSEMBL / HGNC disagree -> use ENSEMBL ID from ENSEMBL
# - ENSEMBL has additional matches -> keep all protein-coding

# I don't remember why I filtered ens_patch for protein_coding but not
# multi_gene

# Manually checked 15/06/22, usually disagreements are because
# HGNC ENSEMBL_ID is obsolete

bm_ids <- dplyr::select(bm, ENSEMBL_ID, HGNC_ID, HGNC_SYMBOL, BIOTYPE) %>%
    unique()

# Genes that differ in HGNC_ID / ENSEMBL_ID combination
bm_ens_diff <- dplyr::anti_join(bm_ids, hgnc_ids,
                                by = c("HGNC_ID", "ENSEMBL_ID"))

# Get the Ensembl IDs that HGNC assigns
hgnc_ens_diff <- hgnc_ids %>%
    AbNames:::filter_by_union(bm_ens_diff %>% dplyr::select(-BIOTYPE))

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
temp <- full_join(temp, multi_gene, relationship = "many-to-many")
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
    dplyr::anti_join(hgnc, by = c("HGNC_SYMBOL", "value")) %>%
    # Only want genes for which there is a HGNC / ENSEMBL /ENTREZ comb in HGNC
    dplyr::semi_join(hgnc %>% dplyr::select(HGNC_ID, ENSEMBL_ID, ENTREZ_ID))

# NCBI and HGNC can differ in how they describe a symbol, e.g. previous
# vs alias. Assume that HGNC annotations are correct, as HGNC is the
# naming consortium.  Remove entries where the value is the same
# (Regardless of what type of symbol it is considered to be)
ncbi_novel <- ncbi %>%
    dplyr::anti_join(hgnc %>% dplyr::select(HGNC_ID, ENSEMBL_ID, value)) %>%
    dplyr::select(-HGNC_NAME)


# UNIPROT IDs can differ.
# Set to HGNC value, as the goal here is not to be comprehensive
hgnc_patch <- hgnc %>%
    dplyr::select(HGNC_SYMBOL, ENSEMBL_ID,ENTREZ_ID, UNIPROT_ID) %>%
    unique()

org_db_novel <- dplyr::rows_update(org_db_novel, hgnc_patch,
                                 by = c("ENSEMBL_ID",
                                        "ENTREZ_ID",
                                        "HGNC_SYMBOL"),
                                 unmatched = "ignore") %>%
    unique()

bm_novel <- dplyr::rows_update(bm_novel, hgnc_patch,
                                   by = c("ENSEMBL_ID",
                                          "ENTREZ_ID",
                                          "HGNC_SYMBOL"),
                                   unmatched = "ignore") %>%
    unique()

# --------------------------------------------------------------------------
# Within HGNC, some aliases differ only in symbol type.
# Prefer symbol over name and current over previous

symbol_order <- c("HGNC_SYMBOL", "HGNC_NAME",
                  "PREVIOUS_SYMBOL", "PREVIOUS_NAME",
                  "ALIAS", "ALIAS_NAME")
hgnc <- hgnc %>%
    mutate(symbol_type = factor(symbol_type, levels = symbol_order)) %>%
    dplyr::group_by(across(-symbol_type)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(symbol_type = as.character(symbol_type))

# --------------------------------------------------------------------------

# Add novel aliases -----

# Add novel aliases from biomaRt and NCBI
# (Note - have previously checked that HGNC_ID / ENSEMBL_ID combinations
# are correct)
gene_aliases <- hgnc %>%
    dplyr::bind_rows(bm_novel, ncbi_novel, org_db_novel) %>%

    # Fill HGNC IDs (missing from org_db)
    dplyr::group_by(ENSEMBL_ID, HGNC_SYMBOL, ENTREZ_ID) %>%
    tidyr::fill(HGNC_ID, UNIPROT_ID, .direction = "updown") %>%

    # Set missing symbol_type to "ALIAS", make "protein-coding" consistent
    dplyr::mutate(symbol_type =
                      ifelse(is.na(symbol_type), "ALIAS", symbol_type),
                  BIOTYPE = ifelse(BIOTYPE == "protein_coding",
                                   "protein-coding", BIOTYPE)) %>%

    # Remove lncRNAs (SIGLEC5 and GOLGA8M, both also have proteins)
    dplyr::filter(! BIOTYPE == "lncRNA") %>%

    # New values may come from more than one source, aggregate
    dplyr::group_by(across(-SOURCE)) %>%
    dplyr::summarise(SOURCE = toString(SOURCE), .groups = "keep") %>%

    # Ensure that aliases are unambiguous
    # (newly added aliases may create new ambiguities)
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_ID)) %>%
    # To check which aliases are ambiguous:
    #dplyr::filter(n_genes > 1 & ! any(symbol_type == "HGNC_SYMBOL"))
    dplyr::filter(n_genes == 1 | symbol_type ==  c("HGNC_SYMBOL")) %>%
    dplyr::select(-n_genes) %>%

    # Aggregate the protein IDs
    dplyr::group_by(HGNC_ID, HGNC_SYMBOL, ENSEMBL_ID, ENTREZ_ID) %>%
    dplyr::mutate(UNIPROT_ID = paste0(sort(unique(UNIPROT_ID)),
                                         collapse = "|"),
                  UNIPROT_ID = dplyr::na_if(UNIPROT_ID, "")) %>%
    ungroup() %>%
    unique()

write_csv(as.data.frame(hgnc), file = "inst/extdata/gene_aliases.csv")
rm(list = setdiff(ls(), c("existing","gene_aliases")))

#gene_aliases <- as.data.frame(hgnc)
#usethis::use_data(gene_aliases, overwrite = TRUE, compress = "bzip2")
