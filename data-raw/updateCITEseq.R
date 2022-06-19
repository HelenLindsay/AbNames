# Add new data to the CITEseq data set and fill missing information if possible
# To do: should we keep any extra columns not in citeseq?
# x: data to add to CITEseq data set
# name: study name
updateCITEseq <- function(x, name = NULL){
    citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
    citeseq <- read.csv(citeseq_fname)

    if (! "Study" %in% colnames(x) & is.null(NULL)){
        stop("Please provide a name for the new data")
    }
    if (! "Antigen" %in% colnames(x)){
        stop("Data to add must contain a column 'Antigen'")
    }

    # Are there already ID columns present?
    id_cols <- c("HGNC_ID", "ENSEMBL_ID", "HGNC_SYMBOL", "Cat_Number")
    id_cols <- intersect(,
                         colnames(x))

    fill_cols <- c()




}








# Bind to CITE-seq

# If all rows have an identifier, group by ID and select majority

# If has catalogue number,
