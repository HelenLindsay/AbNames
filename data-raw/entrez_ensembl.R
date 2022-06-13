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

# ---------------------------------------------------------------------------
# Fill HGNC IDs using org.db if based on symbol and Ensembl


# ---------------------------------------------------------------------------

# Biomart -----

library(tidyverse)
library(biomaRt)
data(hgnc)

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
                  ENTREZ_ID = entrezgene_id) %>%
    dplyr::group_by(HGNC_ID) %>%
    dplyr::mutate(ENTREZ_ID = as.character(ENTREZ_ID),
                  SOURCE = "BIOMART") %>%

    # Only select genes with a HGNC_ID in the hgnc data set
    dplyr::semi_join(hgnc, by = "HGNC_ID") %>%

    # Remove aliases that map to more than one HGNC_SYMBOL
    dplyr::group_by(ALIAS) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_SYMBOL)) %>%
    dplyr::filter(n_genes == 1 | is.na(ALIAS)) %>%
    dplyr::select(-n_genes) %>%
    dplyr::ungroup() %>%

    # Remove entries that don't match HGNC
    # (usually an unofficial symbol, some of these could be saved if necessary)
    dplyr::semi_join(hgnc, by = c("ENSEMBL_ID", "HGNC_SYMBOL"))


# ---------------------------------------------------------------------------
# org.Hs.eg.db genes -----

library(org.Hs.eg.db)

hs <- org.Hs.eg.db

# Select everything that has an Entrez gene ID
org_db <- AnnotationDbi::select(hs,
                 keys = keys(hs),
                 keytype = "ENTREZID",
                 columns = c("SYMBOL", "ALIAS", "ENSEMBL",
                             "GENETYPE", "UNIPROT")) %>%
    dplyr::as_tibble() %>%
    #dplyr::filter(GENETYPE == "protein-coding") %>%
    #dplyr::select(-GENETYPE) %>%
    dplyr::rename(HGNC_SYMBOL = SYMBOL,
                  ENSEMBL_ID = ENSEMBL,
                  ENTREZ_ID = ENTREZID,
                  value = ALIAS,
                  BIOTYPE = gene_biotype) %>%

    # Only keep entries with HGNC symbol (removes e.g. pseudogenes)
    dplyr::semi_join(hgnc, by = "HGNC_SYMBOL") %>%
    dplyr::mutate(SOURCE = "ORGDB") %>%

    # Only keep Ensembl ids that appear in the Biomart table (or missing)
    dplyr::filter(ENSEMBL_ID %in% bm$ENSEMBL_ID | is.na(ENSEMBL_ID)) %>%

    # Remove aliases that map to more than one HGNC_SYMBOL
    dplyr::group_by(value) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_SYMBOL)) %>%
    dplyr::filter(n_genes == 1) %>%
    dplyr::select(-n_genes)

# Fill in the Ensembl ID using hgnc if it is missing, require
# matches between a symbol and (at least) one alias
# Note that semi_join still matches if one value is NA
# (None of the missing genes are in the biomart results)

table(is.na(org_db$ENSEMBL_ID))
# FALSE  TRUE
# 61355   113

hgnc_patch <- hgnc %>%
    dplyr::filter(! is.na(ENSEMBL_ID)) %>%
    dplyr::select(HGNC_SYMBOL, value, ENSEMBL_ID) %>%
    unique()

org_db <- org_db %>%
    dplyr::rows_patch(hgnc_patch,
                      by = c("HGNC_SYMBOL", "value"),
                      unmatched = "ignore")

# None of the remaining NAs can be filled by grouping by HGNC_SYMBOL + ENTREZ_ID
# or by filling from NCBI

table(is.na(org_db$ENSEMBL_ID))
# FALSE  TRUE
# 61400    68


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

# ---------------------------------------------------------------------------
