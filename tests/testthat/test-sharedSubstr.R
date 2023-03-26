test_that(".sharedSubstr correctly groups words", {
    words <- c("fox", "cat", "fox", "in", "box", "fish", "box")
    ids <- c(1, 2, 3, 3, 3, 4, 5)
    expect_equal(.sharedSubstr(words, ids), c(1,2,1,1,1,3,1))
})


test_that("sharedSubstr correctly adds grouping column", {
    df <- data.frame(words = c("fox", "cat", "fox", "in", "box", "fish", "box"),
                     ids = c(1, 2, 3, 3, 3, 4, 5))

    # Should give an error because column names aren't provided
    expect_error(sharedSubstr(df))
    expect_equal(sharedSubstr(df, x = "words", id = "ids"),
                 cbind(df, AB_group = c(1,2,1,1,1,3,1)))
})


# More relevant example df
#df <- data.frame(Antigen = c("CD274 (B7-H1, PD-L1)", "PD-L1", "B7-H1"),
#                 ID = c(1:3))
