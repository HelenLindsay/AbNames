# Add new data to the CITEseq data set and fill missing information if possible
# To do: should we keep any extra columns not in citeseq?
# x: data to add to CITEseq data set
# name: study name
updateCITEseq <- function(x, name = NULL){
    citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
    citeseq <- read.csv(citeseq_fname)

    if (! "Study" %in% colnames(x)){
        if (is.null(name)) {
            stop("Please provide a name for the new data")
        }
        x$Study <- name
    }

    if (! "Antigen" %in% colnames(x)){
        stop("Data to add must contain a column 'Antigen'")
    }

    # Only keep columns already present in CITE-seq collection
    x <- unique(x[, intersect(colnames(x), colnames(citeseq))])

    # CHECK FOR DUPLICATED ANTIGEN NAMES?

    # CHECK FOR NON-ASCI CHARACTERS

    # (DO GREEK SYMBOLS ALSO NEED TO BE REPLACED IN CITESEQ DATA SET?)



    # Are there already ID columns present?
    id_cols <- c("HGNC_ID", "ENSEMBL_ID", "HGNC_SYMBOL", "ENTREZ_ID")
    id_cols <- intersect(id_cols, colnames(x))

    # If ID columns are provided, check that entries are valid
    if (length(id_cols) > 0){

        # ---- to update --------
        data(hgnc)
        # -----------------------

        # Do this as a loop to know which is inconsistent?
        inconsistent <- dplyr::anti_join(x, hgnc, by = id_cols)
        if (nrow(inconsistent) > 0){
            stop("ID rows are inconsistent")
        }
    }

    vendor_info <- c("Cat_Number", "Vendor")


    fill_cols <- c()




}





# Bind to CITE-seq

# If all rows have an identifier, group by ID and select majority

# If has catalogue number,
