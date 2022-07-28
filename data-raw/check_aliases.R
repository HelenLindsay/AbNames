# One HGNC_ID maps to one HGNC_SYMBOL

# One ENTREZ / ENSEMBL ID maps to one HGNC_ID

# One alias maps to one symbol (unless it is a protein complex
# - or multiple isoforms?)

# Coalesce BIOTYPE if this is the only row that differs

# One HGNC_ID to one UNIPROT_ID

# Check that UNIPROT and ENSEMBL IDs are not obsolete
# (any more information in SWISSPROT)

# Coalesce rows = patch with HGNC, group_by + fill missing values
