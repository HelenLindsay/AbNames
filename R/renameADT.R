#'@importFrom methods setGeneric
#'@export
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
#'@export
setMethod("renameADT", as(structure(.Data = c("SingleCellExperiment",
                                              "character"),
                                    names = c("obj", "names"),
                                    package = c("SingleCellExperiment", "")),
                          "signature"),
    function(obj, names, assay = "counts", ...) {
        # to do - add old names to row data

        # names = named vector of new names, names are current names

        # For a SingleCellExperiment, the ADT may either be the main assay
        # or an altExp.
        # May need to replace the row data too
        # sce may have rowPairs, but doesn't appear to inc names
        # toDo: add original names to the rowData
        # check rownames of rowdata

        stopifnot(requireNamespace("SummarizedExperiment"),
                  requireNamespace("SingleCellExperiment"))

        if (assay %in% SingleCellExperiment::altExpNames(obj)){
            # If it's an altExperiment, rename all altExp rownames
            rp_func = altExp
            old_nms <- rownames(SingleCellExperiment::altExp(obj, assay))
            # If there are missing entries, use the old names
            new_nms <- dplyr::coalesce(names[old_nms], old_nms)
            rownames(SingleCellExperiment::altExp(obj, assay)) <- new_nms
            SingleCellExperiment::rowData(
                SingleCellExperiment::altExp(obj, assay)
            ) <- cbind(rowData)

        } else {
            #nms <- .getRownames(obj, assay)

            # Use coalesce in case there are missing values
            new_nms <- dplyr::coalesce(names[rownames(obj)], rownames(obj))
            rownames(obj) <- new_nms
        }
        return(obj)
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
        nms <- rownames(SingleCellExperiment::experiments(obj)[[assay]])
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


# matchToCiteseq -----
# Probably want to standardise name via aliases table.....
# At the moment relies on column names matching
# Rename programmatically?
# Match to the curated citeseq data.set
# nm could be either data.frame or single vector of names
# Eventually, matching just by names may be enough
# To do: check types of matching columns?
# To do: add a column indicating group size?
#
#' Standardise names using the citeseq data set
#'
#'@description Groups antibody names and selects the most frequent name.  By
#'default, the citeseq dataset is grouped by Antigen, Clone, Cat_Number
#' (catalog number) and ALT_ID.
#'@param x A data.frame or tibble containing a column "Antigen" to match to the
#' citeseq data set
#'@param cols (character(n), default NULL) Optional additional columns to use
#'for matching in citeseq data set.  ID columns from the citeseq dataset are
#'used by default.  Antibodies are grouped by a match in any of the ID columns.
#'@param verbose Report which columns are used for matching? (default: TRUE)
#'@param ... Options passed on to getCommonName
#'@export
matchToCiteseq <- function(x, cols = NULL, verbose = TRUE, ...){
    if (! "data.frame" %in% class(x)) {
        stop("x should be a data.frame or tibble")
    }

    cn_x <- colnames(x)

    utils::data("citeseq", envir = environment())

    # Set up names of columns for matching
    keep_cols <- intersect(cols, colnames(citeseq))

    if (! identical(keep_cols, cols)){
        warning("cols should match columns in cite seq data set")
    }

    if (! "Antigen" %in% colnames(x)){
        stop("x must contain a column 'Antigen'")
    }

    # If keep_cols is specified, make sure Antigen is included
    if (! is.null(keep_cols)){
        keep_cols <- unique(c(keep_cols, "Antigen"))
        msg_cols <- toString(keep_cols)

    } else {
        # If columns are not specified, check which of the default columns
        # in getCommon name are present in x
        default_cols <- c("Antigen", "Cat_Number", "Clone", "ALT_ID")
        msg_cols <- toString(intersect(colnames(x), default_cols))
    }

    # Report which columns are used to match new data to cite seq data
    msg <- "Matching new data to citeseq data using columns:\n%s"
    message(sprintf(msg, msg_cols))

    # Add temporary ID column and add new data to citeseq data set
    id <- .tempColName(x, nm = "ID")
    x <- x %>% dplyr::mutate(!!id := "KEEPME")
    x <- dplyr::bind_rows(x, citeseq)

    # If keep_cols is NULL, the default column set is used for matching
    x <- getCommonName(x, cols = keep_cols, ab = "Antigen",
                       fill_col = "Antigen_std", keep = TRUE, ...)

    #######
    # TO DO - ACT ON CHECKS!
    # Check - one Clone should map to one ID
    res <- .checkCiteseq(x, gp = "Cat_Number", id = "HGNC_ID")
    ######

    # Remove ID column, select original columns
    x <- x %>% dplyr::filter(!!sym(id) == "KEEPME") %>%
        dplyr::select(all_of(c("Antigen_std", cn_x, "n_matched")))

    return(x)
}


.checkCiteseq <- function(x, gp = "Cat_Number", id = "HGNC_ID"){
    x %>%
        dplyr::group_by(!!sym(gp)) %>%
        dplyr::filter(dplyr::n_distinct(!!sym(id), na.rm = TRUE) > 1)
}
