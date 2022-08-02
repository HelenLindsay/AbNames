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
setMethod("renameADT", as(structure(.Data = c("SingleCellExperiment",
                                              "character"),
                                    names = c("obj", "names"),
                                    package = c("SingleCellExperiment", "")),
                          "signature"),
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


# Probably want to standardise name via aliases table.....
# At the moment relies on column names matching
# Rename programmatically?
# Match to the curated citeseq data.set
# nm could be either data.frame or single vector of names
# Eventually, matching just by names may be enough
# To do: check types of matching columns?
# To do: add a column indicating group size?
matchToCiteseq <- function(x, cols = NULL){

    # Do do: make citeseq a data set!
    citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
    citeseq <- read.csv(citeseq_fname) %>% unique()

    if (! is.null(cols)){
        keep_cols <- intersect(cols, colnames(citeseq))
        if (! identical(keep_cols, cols)){
            warning("cols should match columns in cite seq data set")
        }
        if (! "data.frame" %in% class(x)){
            stop("If 'cols' argument is provided, x must be a data.frame")
        }
        if (! "Antigen" %in% colnames(x)){
            stop("If 'cols' argument is provided, ",
                 "x must contain a column 'Antigen'")
        }
        cols <- unique(c(keep_cols, "Antigen"))
    }

    if (! "data.frame" %in% class(x)){
        x <- data.frame(Antigen = x)
    }
    if (is.null(cols)) cols <- "Antigen"

    x <- dplyr::bind_rows(x %>% dplyr::mutate(ID = "KEEPME"), citeseq)

    x <- getCommonName(x, cols = cols, ab = "Antigen",
                       new_col = "Antigen_std", keep = TRUE)

    x <- x %>% dplyr::filter(ID == "KEEPME") %>%
        dplyr::select(all_of(c("Antigen", "Antigen_std")))

    # TRIANA HAS DUPLICATES?

}


# getCommonName ----

# default matching columns are Antigen, Cat_Number, Clone, HGNC_SYMBOL)
# cols = columns for grouping
# ab = column for standardising
# new_col = column to add
# ... pass keep = TRUE for keeping grouping columns for debugging
# Be careful of exceptions, e.g. Cat_Number == "custom made"
getCommonName <- function(x, cols = NULL, ab = "Antigen",
                          new_col = "Antigen_std",
                          ignore = list(Cat_Number = "[Cc]ustom"), ...){

    dots <- list(...)

    keep_cols <- c(colnames(x), new_col)

    if (new_col %in% colnames(x)){
        stop(sprintf("Column %s already exists in data.frame", new_col))
    }

    if (is.null(cols)) {
        cols <- c("Antigen", "Cat_Number", "Clone", "HGNC_SYMBOL", "ENSEMBL_ID")
    }

    # Remove sections in brackets, replace Greek symbols
    x <- dplyr::mutate(x, !! new_col := gsub("^[Aa]nti-| \\(.*", "", !!sym(ab)),
                       !! new_col := replaceGreekSyms(!!sym(new_col)))

    # Group by any e.g. catalogue number or exact match to antigen
    tmp_grp <- .tempColName(x, nm = "group")
    x <- group_by_any(x, cols, new_col = tmp_grp, ignore = ignore)

    # Fill with most common value
    x <- fillByGroup(x, group = tmp_grp, method = "all",
                     multiple = "mode", fill = new_col)

    # Problems: Tau (Phospho Thr181) Su and Stephenson?
    # HGNC_SYMBOL WRONG FOR CD158b (KIR2DL2/L3, NKAT2)?  (NKAT2 = only one gene)
    # RRID AB_2810478 matches CD226 and CD98
    # Triana RRID same for CD235a and CD235a - match via combination of RRID
    # and Antigen
    # Triana group 5... - matching because of NA?
    # CD3.1 should not have HGNC symbol PECAM1
    # CD45 MATCHING THROUGH GENE SYMBOL


    if (isTRUE(dots$keep)){
        return(x)
    }

    # Remove temporary columns
    return(dplyr::select(x, all_of(keep_cols)))
}
