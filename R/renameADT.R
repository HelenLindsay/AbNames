#'@importFrom methods setGeneric standardGeneric
setGeneric("renameADT", function(x) {
    standardGeneric("renameADT")
})


# renameADT for signature vector ------
setMethod("renameADT", signature(obj = "character", assay = NULL),
    function(obj, assay){



})



# renameADT for signature SingleCellExperiment ------
#'@importFrom methods setMethod signature
setMethod("renameADT", signature(obj = "SingleCellExperiment",
                                 assay = "character"),
    function(obj, assay, ...) {
        # For a SingleCellExperiment, the ADT may either be the main assay
        # or an altExp.
        # May need to replace the row data too
        # SingleCellExperiment may have rowPairs

        stopifnot(requireNamespace("SummarizedExperiment"),
                  requireNamespace("SingleCellExperiment"))

        if (assay %in% altExpNames(obj)){
            # If it's an altExperiment, rename all altExp rownames
            rp_func = altExp
            nms <- rownames(altExp(obj, assay))
        } else {
            nms <- .getRownames(obj, assay)
        }


})


# renameADT for signature MultiAssayExperiment ------
setMethod("renameADT", signature(obj = "MultiAssayExperiment",
                                 assay = "character"),
    function(obj, assay, ...) {
        stopifnot(requireNamespace("SummarizedExperiment"),
                  requireNamespace("MultiAssayExperiment"))
        nms <- rownames(experiments(obj)[[assay]])

        #nms <- .getRownames(obj, assay)

})


# renameADT for signature Seurat ------
setMethod("renameADT", signature(obj = "Seurat",
                                 assay = "character"),
    function(obj, assay) {

})


# Helper functions -----
.getRownames <- function(obj, assay, ...){
    if (! assay %in% names(assays(obj))){
        stop(sprintf("Assay %s not found", assay))
    }
    return(SummarizedExperiment::assays(obj)[[assay]])
}
