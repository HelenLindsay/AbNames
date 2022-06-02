# Notes -----

# Ensembl is based on Refseq, Uniprot and Havana

# One Ensembl may map to multiple EntrezGene if there are multiple good, but no
# perfect, matches.

# Neither Ensembl nor REFSEQ protein/transcript identifiers appear to be
# matched to isoform names, e.g. CD45 -> CD45RO

# Fetch org.Hs.eg.db genes -----

library(org.Hs.eg.db)

hs <- org.Hs.eg.db

# Select everything that has an ensembl protein ID
org_db <- select(hs,
              keys = mappedkeys(org.Hs.egENSEMBLPROT),
              keytype = "ENTREZID",
              columns = c("SYMBOL", "ALIAS", "ENSEMBL")) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(HGNC_SYMBOL = SYMBOL,
                  ENSEMBL_ID = ENSEMBL,
                  ENTREZ_ID = ENTREZID)

