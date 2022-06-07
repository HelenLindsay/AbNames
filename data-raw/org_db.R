# Notes -----

# Ensembl is based on Refseq, Uniprot and Havana

# One Ensembl may map to multiple EntrezGene if there are multiple good, but no
# perfect, matches.

# Neither Ensembl nor REFSEQ protein/transcript identifiers appear to be
# matched to isoform names, e.g. CD45 -> CD45RO

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
    dplyr::mutate(n_entrez = n_distinct(ENTREZ_ID),
                  n_ensembl = n_distinct(ENSEMBL_ID),
                  ENTREZ_ID = as.character(ENTREZ_ID)) %>%
    # Only select genes with a HGNC_ID in the hgnc data set
    dplyr::semi_join(hgnc, by = "HGNC_ID")



# Fetch org.Hs.eg.db genes -----

library(org.Hs.eg.db)

hs <- org.Hs.eg.db

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
    # Only keep Ensembl ids that appear in the Biomart table
    dplyr::filter(ENSEMBL_ID %in% bm$ENSEMBL_ID)

# Does org_db have anything not in NCBI genes (and vice versa)? ----

org_db_long <- org_db %>%
    dplyr::rename(NCBI_ID = ENTREZ_ID) %>%
    tidyr::pivot_longer(cols = c(HGNC_SYMBOL, ALIAS),
                        names_to = "symbol_type") %>%
    dplyr::select(-symbol_type) %>%
    unique()

ncbi_genes_long <- ncbi_genes %>%
    dplyr::mutate(NCBI_ID = as.character(NCBI_ID)) %>%
    tidyr::pivot_longer(cols = c(HGNC_SYMBOL, value),
                        names_to = "symbol_type2") %>%
    dplyr::select(-symbol_type, -symbol_type2) %>%
    unique()

# Yes, org_db has aliases that are not in ncbi_genes
org_db_not_ncbi <- org_db_long %>%
    dplyr::anti_join(ncbi_genes_long,
                     by = c("NCBI_ID", "ENSEMBL_ID", "value"))


# NCBI genes includes CD antigen descriptions, e.g. "CD138 antigen"
ncbi_not_orgdb <- ncbi_genes_long %>%
    dplyr::anti_join(org_db_long,
                     by = c("NCBI_ID", "ENSEMBL_ID", "value"))


# Check if ncbi_genes and org_db agree on mapping between
# HGNC, ENTREZ and ENSEMBL -----






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

