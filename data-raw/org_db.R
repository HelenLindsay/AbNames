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
                  ENTREZ_ID = ENTREZID)

