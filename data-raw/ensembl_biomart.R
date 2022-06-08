library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

chrs <- c(1:22, "MT", "X", "Y")

bm <- getBM(mart = ensembl,

            filters = c("ensembl_gene_id",
                        "chromosome_name",
                        "with_hgnc",
                        "biotype"),

            values = list("ENSG00000081237",
                          chrs,
                          TRUE,
                          "protein_coding"),

            attributes = c("ensembl_gene_id",
                           "description",
                           "external_gene_name",
                           "external_gene_source",
                           "gene_biotype",
                           "hgnc_id",
                           "hgnc_symbol")
            )
