test_that("renameADT for SingleCellExperiment works", {
    # SingleCellExperiment with ADT as main experiment

    u <- matrix(rpois(2000, 5), ncol=100)
    v <- matrix(rpois(20, 5), ncol=100)

    # Some non-standard names, including duplicates
    adt_nms <- c("KLRG1 (MAFA)", "KLRG1", "CD3 (CD3E)", "HLA.A.B.C",
                 "HLA-A/B/C", "NKAT2", "CD66a.c.e", "CD66a_c_e",
                 "CD11a/CD18 (LFA-1)", "")

    # NKAT2 = CD158b
    test_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts=u),
        altExps = list(adt = SingleCellExperiment::SingleCellExperiment(v)))

    # SingleCellExperiment with ADT as altExp
})
