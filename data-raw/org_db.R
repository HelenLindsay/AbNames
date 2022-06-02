library(org.Hs.eg.db)

hs <- org.Hs.eg.db

# Select everything that has an ensembl protein ID
org_db <- select(hs,
              keys = mappedkeys(org.Hs.egENSEMBLPROT),
              keytype = "ENTREZID",
              columns = c("SYMBOL", "ALIAS",
                          "ENSEMBL", "ENSEMBLPROT", "REFSEQ")) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(HGNC_SYMBOL = SYMBOL,
                  ENSEMBL_ID = ENSEMBL,
                  ENTREZ_ID = ENTREZID)

