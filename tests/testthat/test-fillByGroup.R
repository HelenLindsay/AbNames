test_that("fillByGroup stops with multiple values", {
    df <- data.frame(A = c(rep("A", 3), NA, rep("B", 4), rep("C", 3), NA),
                     B = c(rep("A", 2), NA, "A", rep("B", 4), rep("C", 3), NA),
                     C = c(NA, 1, 1, 1, 2, 3, NA, 2, 4, 4, NA, 5))

})
