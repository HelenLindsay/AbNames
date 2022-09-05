test_that("sharedSubstr correctly groups words", {
    words <- c("fox", "cat", "fox", "in", "box", "fish", "box")
    ids <- c(1, 2, 3, 3, 3, 4, 5)
    expect_equal(.sharedSubstr(words, ids), c(1,2,1,1,1,3,1))
})
