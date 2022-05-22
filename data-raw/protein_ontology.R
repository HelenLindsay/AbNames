# Download the protein ontology table.
# Here we need the synonyms but don't need the ontology relationships as these
# are available in the package ontoProc

# Test cases CD77 and CD66

library(ontologyIndex)

fn <- "https://proconsortium.org/download/current/pro_reasoned.obo"
f <- tempfile()
download.file(fn, destfile = f)

po <- ontologyIndex::get_ontology(f,
                                  propagate_relationships =
                                    c("is_a", "intersection_of","relationship"),
                                  extract_tags = "everything")
id <- po$id[! grepl("Araport|CGNC|EcoGene|EnsemblBacteria|FlyBase", po$id)]

# NCBITaxon:9606

# BFO = biological process ontology of GO
# HGNC
# NCBIGene
# NCBITaxon
# Ensembl
