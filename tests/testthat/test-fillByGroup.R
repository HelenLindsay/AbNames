df <- data.frame(A = c(rep("A", 3), NA, rep("B", 4), rep("C", 3), NA),
                 B = c(rep("A", 2), NA, "A", rep("B", 4), rep("C", 3), NA),
                 C = c(NA, 1, 1, 1, 2, 3, NA, 2, 4, 4, NA, 5))


test_that("fillByGroup with option=stop stops with multiple values", {
    expect_error(fillByGroup(df, multiple = "stop"))
})

test_that("fillByGroup with option=majority fills multiple values", {
    exp_res <- data.frame(A = c(rep("A", 2), rep("B", 4),
                                rep("C", 3), "A", NA),
                          B = c(rep("A", 2), rep("B", 4),
                                rep("C", 3), rep(NA, 2)),
                     C = c(1, 1, 2, 3, 2, 2, 4, 4, 4, 1, 5))

   expect_equal(fillByGroup(df, multiple = "majority"), exp_res)

})
