# Notes -----

# Ensembl is based on Refseq, Uniprot and Havana

# One Ensembl may map to multiple EntrezGene if there are multiple good, but no
# perfect, matches.

# Neither Ensembl nor REFSEQ protein/transcript identifiers appear to be
# matched to isoform names, e.g. CD45 -> CD45RO

# org.Hs.egENSEMBLPROT doesn't include all protein-coding genes, e.g.
# ENSG00000125378 BMP4 is not matched to an ENSEMBLPROT ID

# Some Ensembl genes do not have a match in org.db, e.g. Entrez 9103 = FCGR2C
# ENSG00000244682 - (From Ensembl website: does not contain any transcripts for
# which we have selected identical models in RefSeq).

# Some Ensembl genes map to multiple Entrez genes, e.g. ENSG00000276070 to
# 9560 / 388372

# Some Ensembl genes map to multiple HGNC_IDs, e.g. ENSG00000276085 to
# HGNC:10628 / HGNC:30554

# Note that fetching chromosomes only misses haplotypes, which include many
# immune genes

# ---------------------------------------------------------------------------
existing <- ls()
hgnc <- read_delim(file = "inst/extdata/hgnc.csv")

# Biomart -----

# Note: adding uniprotswissprot as attribute acts like a filter,
# e.g. "ENSG00000105501" is not present in output

chrs <- c(1:22, "X", "Y", "MT")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

bm <- getBM(mart = ensembl,
          filters = c("chromosome_name", "with_hgnc"),
          values = list(chrs, TRUE),
          attributes = c("ensembl_gene_id", "external_gene_name",
                         "external_synonym", "hgnc_id", "entrezgene_id",
                         "gene_biotype")) %>%
    tibble::as_tibble() %>%
    dplyr::rename(ENSEMBL_ID = ensembl_gene_id,
                  HGNC_SYMBOL = external_gene_name,
                  ALIAS = external_synonym,
                  HGNC_ID = hgnc_id,
                  ENTREZ_ID = entrezgene_id,
                  BIOTYPE = gene_biotype) %>%
    dplyr::mutate(ENTREZ_ID = as.character(ENTREZ_ID)) %>%
    # Remove pseudogenes
    dplyr::filter(! grepl("pseudo", BIOTYPE))


# Add the UNIPROT IDs (SWISSPROT ONLY) ----

ens_to_sp <- getBM(mart = ensembl,
            filters = c("chromosome_name", "with_hgnc"),
            values = list(chrs, TRUE),
            attributes = c("ensembl_gene_id", "hgnc_id", "external_synonym",
                           "external_gene_name","uniprotswissprot")) %>%
    dplyr::rename(ENSEMBL_ID = ensembl_gene_id,
                  UNIPROT_ID = uniprotswissprot,
                  HGNC_SYMBOL = external_gene_name,
                  HGNC_ID = hgnc_id,
                  ALIAS = external_synonym) %>%
    dplyr::mutate(UNIPROT_ID = na_if(UNIPROT_ID, "")) %>%
    dplyr::filter(! is.na(UNIPROT_ID))

if (! all(ens_to_sp$ENSEMBL_ID %in% bm$ENSEMBL_ID)){
    warning("Uniprot table contains genes missing from bm table")
}

# Rows will be gained as one gene may have multiple UNIPROT IDs
bm <- bm %>%
    dplyr::full_join(ens_to_sp) %>%
    # Only select genes with a HGNC_ID in the hgnc data set
    dplyr::semi_join(hgnc, by = "HGNC_ID")


# Fix the Biomart entries that have the wrong official symbol ----
hgnc_ids <- hgnc %>%
    dplyr::select(ENSEMBL_ID, HGNC_ID, HGNC_SYMBOL) %>%
    unique()

bm_patch <- bm %>%
    # Ensembl ID and HGNC ID are shared
    dplyr::semi_join(hgnc, by = c("ENSEMBL_ID", "HGNC_ID")) %>%
    # But symbol is incorrect
    dplyr::anti_join(hgnc, by = c("ENSEMBL_ID", "HGNC_SYMBOL", "HGNC_ID"))
# length(unique(bm_patch$ENSEMBL_ID)) 66 genes

# Remove these rows from biomart
bm <- bm %>% anti_join(bm_patch)

# If the incorrect symbol is a known alias in HGNC, fix and add back in
bm_patch <- bm_patch %>%
    # If the correct symbol is a known HGNC alias ...
    dplyr::semi_join(hgnc, by = c("ENSEMBL_ID", "HGNC_SYMBOL" = "value")) %>%
    # Make the incorrect symbol an ALIAS
    dplyr::group_by(ENSEMBL_ID) %>%
    tidyr::pivot_longer(cols = c(HGNC_SYMBOL, ALIAS), values_to = "ALIAS",
                        names_to = NULL) %>%
    unique() %>%
    # Add the HGNC_SYMBOL from hgnc table
    dplyr::left_join(hgnc_ids)

# Add patched rows back into Biomart
bm <- bm %>% dplyr::bind_rows(bm_patch)

# --------

# Remove aliases that map to more than one HGNC_SYMBOL ----
bm <- bm %>%
    dplyr::mutate(ALIAS = na_if(ALIAS, "")) %>%
    dplyr::group_by(ALIAS) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_SYMBOL)) %>%
    dplyr::filter(n_genes == 1 | is.na(ALIAS)) %>%
    dplyr::select(-n_genes) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SOURCE = "BIOMART")

# ---------------------------------------------------------------------------
# org.Hs.eg.db genes -----

hs <- org.Hs.eg.db

# Select everything that has an Entrez gene ID
org_db <- AnnotationDbi::select(hs,
                 keys = keys(hs),
                 keytype = "ENTREZID",
                 columns = c("SYMBOL", "ALIAS", "ENSEMBL",
                             "GENETYPE", "UNIPROT")) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(HGNC_SYMBOL = SYMBOL,
                  ENSEMBL_ID = ENSEMBL,
                  ENTREZ_ID = ENTREZID,
                  value = ALIAS,
                  BIOTYPE = GENETYPE,
                  UNIPROT_ID = UNIPROT) %>%

    # Only keep entries with HGNC symbol (removes e.g. pseudogenes)
    dplyr::semi_join(hgnc, by = "HGNC_SYMBOL") %>%
    dplyr::mutate(SOURCE = "ORGDB") %>%

    # Only keep Ensembl ids that appear in the Biomart table (or missing),
    # remove pseudogenes by biotype
    dplyr::filter(ENSEMBL_ID %in% bm$ENSEMBL_ID | is.na(ENSEMBL_ID),
                  ! grepl("pseudo|unknown|ncRNA", BIOTYPE)) %>%

    # Remove aliases that map to more than one HGNC_SYMBOL
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_SYMBOL)) %>%
    dplyr::filter(n_genes == 1) %>%
    dplyr::select(-n_genes)

# Fill in the Ensembl ID using hgnc if it is missing, require ---
# matches between a symbol and (at least) one alias
# Note that semi_join still matches if one value is NA
# (None of the missing genes are in the biomart results)

table(is.na(org_db$ENSEMBL_ID))
# FALSE  TRUE
# 117602    206


# Patch Ensembl ID by symbol / Entrez ID
hgnc_patch <- hgnc %>%
    dplyr::filter(! is.na(ENSEMBL_ID)) %>%
    dplyr::select(HGNC_SYMBOL, ENSEMBL_ID, ENTREZ_ID) %>%
    unique()

org_db <- org_db %>%
    dplyr::rows_patch(hgnc_patch,
                      by = c("HGNC_SYMBOL", "ENTREZ_ID"),
                      unmatched = "ignore")


#table(is.na(org_db$ENSEMBL_ID))
# FALSE   TRUE
# 117706    106

