test_that("renameADT for SingleCellExperiment works", {
    # SingleCellExperiment with ADT as main experiment

    # gene expression
    u <- matrix(rpois(5000, 5), ncol=100)
    rownames(u) <- sprintf("GENE%s", 1:50)

    # adt expression
    v <- matrix(rpois(2000, 5), ncol=100)
    rownames(v) <- sprintf("CD%s", 1:20)
    new_adt_nms <- sprintf("ADT%s", 1:20)

    # Another dummy altExp, to check that names aren't changed
    w <- matrix(rpois(2000, 5), ncol=100)
    rownames(w) <- sprintf("DONT_CHANGE_%s", 1:20)

    adt_sce <- SingleCellExperiment::SingleCellExperiment(v)
    dummy_sce <- SingleCellExperiment::SingleCellExperiment(w)
    test_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts=u), altExps = list(adt = adt_sce))

    # No assay name, only one assay
    res1 <- renameADT(adt_sce, new_adt_nms)
    expect_equal(rownames(SummarizedExperiment::assay(res1)), new_adt_nms)

    # SingleCellExperiment with ADT as altExp, non-default assay name
    res2 <- renameADT(test_sce, new_adt_nms, assay = "adt")
    # Rownames of genes shouldn't change
    expect_equal(rownames(SummarizedExperiment::assay(res2)), rownames(u))
    # Rownames of ADT should be updated
    expect_equal(rownames(SingleCellExperiment::altExp(res2, "adt")),
                 new_adt_nms)

    # Add a dummy experiment same size as ADT
    SingleCellExperiment::altExp(test_sce, "dummy") <- dummy_sce

    # Expect an error if we don't specify the assay
    expect_error(renameADT(test_sce, new_adt_nms))

    res3 <- renameADT(test_sce, new_adt_nms, assay = "adt")
    # Rownames of genes shouldn't change
    expect_equal(rownames(SummarizedExperiment::assay(res3)), rownames(u))
    # Rownames of ADT should be updated
    expect_equal(rownames(SingleCellExperiment::altExp(res3, "adt")),
                 new_adt_nms)
    # Rownames of dummy altExp shouldn't change
    expect_equal(rownames(SingleCellExperiment::altExp(res3, "dummy")),
                 rownames(w))

    # To do: check that original names are recorded
})
