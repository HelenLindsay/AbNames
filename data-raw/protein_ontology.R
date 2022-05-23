# Download the protein ontology table.
# Here we need the synonyms but don't need the ontology relationships as these
# are available in the package ontoProc

# Test cases CD77 and CD66

library(ontologyIndex)
library(tidyverse)

fn <- "https://proconsortium.org/download/current/pro_reasoned.obo"
f <- tempfile()
download.file(fn, destfile = f)

po <- ontologyIndex::get_ontology(f,
                                  #propagate_relationships =
                                  #  c("is_a", "intersection_of","relationship"),
                                  extract_tags = "everything")



# Just get the entries with a synonym
has_synonym <- lengths(po$synonym) > 0
po_syns <- do.call(cbind, po)[has_synonym,]
po_syns <- tibble::as_tibble(po_syns) %>%
    tidyr::unnest(cols = c(id, name, obsolete, only_in_taxon),
                  keep_empty = TRUE) %>%
    # Remove proteins that aren't in human
    dplyr::filter(is.na(only_in_taxon) | only_in_taxon == "NCBITaxon:9606") %>%
    dplyr::select(id, name, synonym, obsolete, xref) %>%
    tidyr::unnest(cols = synonym, keep_empty = TRUE) %>%
    dplyr::mutate(synonym_id = gsub(".* \\[(.*)\\]$", "\\1", synonym),
                  synonym_type = gsub('\\".*\\" (.*) \\[.*$', "\\1", synonym),
                  synonym = gsub('\\"(.*)\\".*', "\\1", synonym),
                  uniprot_id = .gsubNA("(UniProt[A-z0-9:]+),.*", "\\1",
                                       synonym)) %>%
    tidyr::unnest(xref, keep_empty = TRUE)


# MOD: ID for protein modification

# Only interested in entries that have an id?



    # Remove columns with no information
    dplyr::select(-`auto-generated-by`, -`data-version`, -date,
                  -`default-namespace`, -domain, -`format-version`, -idspace,
                  -is_symmetric, -ontology, -property_value, -range,
                  -`saved-by`, -synonymtypedef)


# Keep:
# c(name, id, synonym, def, alt_id, has_component, is_a, has_part,
# intersection_of, involved_in, replaced_by, union_of)

# Is "EXACT" in synonym always the same as xref?

# Check which columns have no information:
#lapply(po_syns, function(x) unique(lengths(x)))


# BFO = biological process ontology of GO
# HGNC
# NCBIGene
# NCBITaxon
# Ensembl
