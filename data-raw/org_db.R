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

# Biomart -----

library(biomaRt)
data(hgnc)

chrs <- c(1:22, "X", "Y", "MT")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

bm <- getBM(mart = ensembl,
            filters = c("chromosome_name", "with_hgnc", "biotype"),
            values = list(chrs, TRUE, "protein_coding"),
            attributes = c("ensembl_gene_id", "external_gene_name",
              "external_synonym", "hgnc_id", "entrezgene_id")) %>%
    tibble::as_tibble() %>%
    dplyr::rename(ENSEMBL_ID = ensembl_gene_id,
                  HGNC_SYMBOL = external_gene_name,
                  ALIAS = external_synonym,
                  HGNC_ID = hgnc_id,
                  ENTREZ_ID = entrezgene_id) %>%
    dplyr::group_by(HGNC_ID) %>%
    dplyr::mutate(ENTREZ_ID = as.character(ENTREZ_ID),
                  source = "BIOMART") %>%

    # Only select genes with a HGNC_ID in the hgnc data set
    dplyr::semi_join(hgnc, by = "HGNC_ID") %>%

    # Remove aliases that map to more than one HGNC_SYMBOL
    dplyr::group_by(ALIAS) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_SYMBOL)) %>%
    dplyr::filter(n_genes == 1) %>%
    dplyr::select(-n_genes) %>%

    # Remove entries that don't match HGNC (usually an unofficial symbol)
    dplyr::semi_join(hgnc, by = c("ENSEMBL_ID", "HGNC_SYMBOL"))


# Select the aliases that aren't already present in hgnc
bm_novel <- bm %>%
    dplyr::mutate(symbol_type = "ALIAS") %>%
    dplyr::rename(value = ALIAS) %>%
    dplyr::filter(! value == "") %>%
    dplyr::anti_join(hgnc, by = c("HGNC_SYMBOL", "value"))



# Fetch org.Hs.eg.db genes -----

library(org.Hs.eg.db)

hs <- org.Hs.eg.db

# Select everything that has an Entrez gene ID
org_db <- select(hs,
                 keys = keys(hs),
                 keytype = "ENTREZID",
                 columns = c("SYMBOL", "ALIAS", "ENSEMBL")) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(HGNC_SYMBOL = SYMBOL,
                  ENSEMBL_ID = ENSEMBL,
                  ENTREZ_ID = ENTREZID) %>%

    # Only keep entries with HGNC symbol (removes e.g. pseudogenes)
    dplyr::semi_join(hgnc, by = "HGNC_SYMBOL") %>%
    dplyr::mutate(source = "ORGDB") %>%

    # Only keep Ensembl ids that appear in the Biomart table (or missing)
    dplyr::filter(ENSEMBL_ID %in% bm$ENSEMBL_ID | is.na(ENSEMBL_ID))

    # Remove aliases that map to more than one HGNC_SYMBOL
    dplyr::group_by(ALIAS) %>%
    dplyr::mutate(n_genes = n_distinct(HGNC_SYMBOL)) %>%
    dplyr::filter(n_genes == 1) %>%
    dplyr::select(-n_genes)

# Select the aliases that aren't already present in hgnc -----
org_db_novel <- org_db %>%
    dplyr::mutate(symbol_type = "ALIAS") %>%
    dplyr::rename(value = ALIAS) %>%
    dplyr::anti_join(hgnc, by = c("HGNC_SYMBOL", "value"))


# Select the ncbi aliases that aren't already present in hgnc ----
ncbi_novel <- ncbi_genes %>%

    dplyr::mutate(symbol_type = "ALIAS") %>%
    dplyr::rename(value = ALIAS) %>%
    dplyr::anti_join(hgnc, by = c("HGNC_SYMBOL", "value"))










# Check if ncbi_genes and org_db agree on mapping between
# HGNC, ENTREZ and ENSEMBL -----

ncbi_id_map <- ncbi_genes %>%
    dplyr::select(NCBI_ID, ENSEMBL_ID) %>%
    unique()

org_db_id_map <-  org_db %>%
    dplyr::select(ENTREZ_ID, ENSEMBL_ID) %>%
    unique()

ncbi_diff <- anti_join(ncbi_id_map, org_db_id_map)
org_db_diff <- anti_join(org_db_id_map, ncbi_id_map)



# Joining Entrez and Ensembl -----

# Usually one ENSEMBL to one HGNC (occasional exceptions), e.g. TCBE
# TCBE appears to be multiple copies that produce the same protein
# ENTREZ can be several to one HGNC

# Within org_db, one ENTREZ_ID can map to several ENSEMBL_IDs
# e.g. APOBEC3A_B / APOBEC3A BiomaRt only maps to one of the aliases



# Sometimes biomaRt ENTREZ_IDs are NA.
# Fill these if the ENSEMBL_ID and HGNC_SYMBOL agree

org_db_patch <- org_db %>%
    dplyr::select(ENSEMBL_ID, HGNC_SYMBOL, ENTREZ_ID) %>%
    dplyr::filter(! is.na(ENSEMBL_ID)) %>%
    unique()

bm <- bm %>%
    dplyr::rows_patch(org_db_patch,
                       by = c("ENSEMBL_ID", "HGNC_SYMBOL"),
                      unmatched = "ignore")






union_join <- function(df, rows){
    x <- unlist(df[rows, ], use.names = FALSE)
    df %>% dplyr::filter(if_any(.cols = everything(), ~.x %in% x))
}




# There are cases where Biomart and org_db map the same
# ENSEMBL_ID / HGNC_SYMBOL to different ENTREZ_IDs,
# e.g. ENSG00000143702/CEP170 to 645455 (biomaRt) or 9859 (org_db)
# In this case, biomaRt mapping is to a pseudogene according to NCBI


# Just keep rows from org_db and Biomart where Entrez / Ensembl agree

bm_minus_entrez <- dplyr::anti_join(bm, org_db,
                                    by = c("ENSEMBL_ID", "ENTREZ_ID"))

ens_entrez <- dplyr::semi_join(bm, org_db,
                               by = c("ENSEMBL_ID", "ENTREZ_ID"))


# Two problems: one is mapping the IDs, the other is collecting the aliases
