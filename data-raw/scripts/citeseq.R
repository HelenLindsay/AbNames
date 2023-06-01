# Create citeseq data set ----
# The citeseq data set is a combined table of human antibodies
# used in ~50 publically available studies.
# The steps involved in creating the merged clone table are documented in
# github.com/HelenLindsay/CITEseq_antibody_data

citeseq <- readr::read_delim("../raw_data/merged_adt_clones.tsv")

citeseq <- as.data.frame(citeseq)
usethis::use_data(citeseq, overwrite = TRUE, compress = "bzip2")
