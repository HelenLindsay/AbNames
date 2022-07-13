#'@importFrom methods setGeneric
setGeneric("renameADT", signature = c("obj", "names"),
           function(obj, names, ...) {
    standardGeneric("renameADT")
})


# renameADT for signature vector ------
setMethod("renameADT", signature(obj = "character", names = "character"),
    function(obj, names){

})



# renameADT for signature SingleCellExperiment ------
#'@importFrom methods setMethod signature
setMethod("renameADT", signature(obj = "SingleCellExperiment",
                                 names = "character"),
    function(obj, names, assay = "counts", ...) {

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
# https://stackoverflow.com/questions/57380044/
# how-to-document-s4-methods-that-rely-on-classes-from-external-packages
setMethod("renameADT", as(structure(.Data = c("MultiAssayExperiment",
                                              "character"),
                                     names = c("obj", "names"),
                                     package = c("MultiAssayExperiment", "")),
                          "signature"),
    function(obj, names, ...) {
        stopifnot(requireNamespace("SummarizedExperiment"),
                  requireNamespace("MultiAssayExperiment"))
        nms <- rownames(experiments(obj)[[assay]])
        #nms <- .getRownames(obj, assay)
})


# # renameADT for signature Seurat ------
# #'@importFrom methods slotNames slot
# setMethod("renameADT", signature(obj = "Seurat",
#                                  assay = "character"),
#     function(obj, assay) {
#         stopifnot(requireNamespace("Seurat"),
#                   requireNamespace("SeuratObject"))
#
#         if (! assay %in% Assays(obj)){
#             stop(sprintf("Assay %s not found.\n  Available assays: %s\n",
#                          assay, toString(SeuratObject::Assays(obj))))
#         }
#
#         # To do: add Seurat, SeuratObject to suggests
#         # Names might also have to be changed in "reductions"?
#         # (not in Triana_2022, maybe in others)
#
#         # Can rownames be updated just with rownames?
#         # Or convert to SingleCellExperiment and back?
#
#
#         # Find slots in the Assay object that have to be changed
#         sr_assay <- obj@assays[[assay]]
#         assay_slots <- slotNames(sr_assay)
#         assay_has_rn <- vapply(assay_slots, function(nm){
#             ! is.null(rownames(slot(sr_assay, nm)))
#             }, logical(1))
#
#         # https://github.com/satijalab/seurat/issues/1049
#         # https://github.com/mojaveazure/seurat-object/issues/3
#         #obj@assays$RNA@counts@Dimnames[[1]] <- mgi$genes
#         #obj@assays$RNA@data@Dimnames[[1]] <- mgi$genes
#         #obj@assays$RNA@meta.features <- mgi
#
#         # modify counts, data and scale.data slots for every assay
#         # you have in your Seurat object. In addition, you will also
#         # need to modify the meta.data slots for both the object and
#         # each assay object contained within the Seurat object.
#         # Finally, you'll need to rerun any dimensional reduction
#         # and graph-building (eg. FindNeighbors)
#
#         # RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
#         #print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
#         #RNA <- obj@assays$RNA
#
#         #if (nrow(RNA) == length(newnames)) {
#         #    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]    <- newnames
#         #    if (length(RNA@data)) RNA@data@Dimnames[[1]]        <- newnames
#         #    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]
#         # <- newnames
#         #} else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
#         #obj@assays$RNA <- RNA
#         #return(obj)
#     #}
#
# })


# Helper functions -----
.getRownames <- function(obj, assay, ...){
    if (! assay %in% names(SummarizedExperiment::assays(obj))){
        stop(sprintf("Assay %s not found", assay))
    }
    return(SummarizedExperiment::assays(obj)[[assay]])
}
