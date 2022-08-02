test_that("upperSquish correctly ignores groups of numbers", {
    test <- c("CD3.1", "CD3-A", "Nectin-2", "IL-2Rb", "CD8a", "LFA.3")
    exp <- c(NA, "CD3A", "NECTIN2", "IL2RB", "CD8A", "LFA3")
    expect_equal(upperSquish(test), exp)
})
