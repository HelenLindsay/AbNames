test_that("splitMerge gives the same answer with quoted or unquoted filter", {
    df <- data.frame(A = c(1, 1, 1, 2, NA, 3, NA),
                     B = c(4, 4, NA, 5, 6, 7, NA), C = 8:14)

    f <- function(df, x) df %>% dplyr::mutate(D = {{ x }} > 10)

    res_unquoted <- splitMerge(df, complete.cases(A, B), f, x=C)
    res_quoted <- splitMerge(df, "complete.cases(A, B)", f, x=C)
    expect_equal(res_unquoted, res_quoted)
    # Check the result is actually correct:
    expect_equal(res_quoted$D, c(rep(FALSE, 2), rep(TRUE, 2), rep(NA, 3)))
})


test_that("splitMerge works with multiple filters", {
    df <- data.frame(A=c(1, 1, 1, 2, NA, 3, NA),
                     B=c(4, 4, NA, 5, 6, 7, NA), C=8:14)
    f <- function(df, x) dplyr::mutate(df, D = {{ x }} + 1)

    res <- splitMerge(df, "A < 3 & B == 4", f, x=C)
    expect_equal(res, dplyr::mutate(df, D = c(9, 10, rep(NA, 5))))
})
