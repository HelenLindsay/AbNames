test_that(".tempColName returns a column name that doesn't already exist", {
    expect_equal(.tempColName(data.frame(x = 1)), "TEMP")
    expect_equal(.tempColName(data.frame(TEMP = 1)), "TEMP1")
    expect_equal(.tempColName(data.frame(x = 1), n = 3),
                 c("TEMP", "TEMP1", "TEMP2"))
    expect_equal(.tempColName(data.frame(TEMP = 1, TEMP2 = 2), n = 2),
                 c("TEMP1", "TEMP3"))
})


test_that(".wordGrep handles pattern at all positions", {
    expect_equal(.wordGrep("CD3", "CD33"), FALSE)
    expect_equal(.wordGrep("CD3", "in the CD3 middle"), TRUE)
    expect_equal(.wordGrep("CD3", "sentence about CD3."), TRUE)
    expect_equal(.wordGrep("CD3", "two CD3. Sentences."), TRUE)
    expect_equal(.wordGrep("CD3", "CD3"), TRUE)
    expect_equal(.wordGrep("CD3", "CD3.A"), FALSE)
    expect_equal(.wordGrep("CD3", "Sentence with CD3.A"), FALSE)
})


test_that(".stopIfColExists is vectorised", {
    df <- setNames(data.frame(as.list(1:3)), LETTERS[1:3])
    expect_error(.stopIfColExists(df, "A"),
                 regexp = "^Column(s) A already.*|Please.*")
    expect_error(.stopIfColExists(df, c("A", "B")),
                 regexp = "Column(s) A, B already|Please.*")
    expect_error(.stopIfColExists(df, c("A", "B", "E")),
                 regexp = "Column(s) A, B already|Please.*")
    expect_equal(.stopIfColExists(df, c("E")), NULL)
})


test_that("union_join works as expected", {
    data(diamonds)
    diamonds <- diamonds[1:20,] %>%
        dplyr::mutate(across(c("cut", "color", "clarity"), as.character))
    res1 <- union_join(diamonds, data.frame(cut = "Fair", color = "I"))
    exp_res <- diamonds %>% dplyr::filter(color == "I" | cut == "Fair")
    expect_equal(res1, exp_res)

    # Value 4.05 appears in columns x and y.
    # union_join should respect column
    res2 <- union_join(diamonds, data.frame(x = 4.05))
    exp_res2 <- diamonds %>% dplyr::filter(x == 4.05)
    expect_equal(res2, exp_res2)

    # If one data.frame df and a row range is supplied, union_join should
    # subset from df
    res3 <- union_join(diamonds, rows = 10)
    exp_res3 <- diamonds %>%
        dplyr::filter(carat == 0.23 | cut == "Very Good" | color == "H" |
                          clarity == "VS1" | table == 61 )
    expect_equal(res3, exp_res3)

    # If two data.frames and a row.range are supplied, union_join should
    # subset from second data.frame
    res4 <- union_join(diamonds, diamonds[, c("cut", "color")], rows = 9:10)
    res5 <- union_join(diamonds, diamonds[9:10, c("cut", "color")])
    expect_equal(res4, res5)

    # Using "by" to subset should be equivalent to explicitly subsetting df2
    res6 <- union_join(diamonds, diamonds[9:10,], by = c("cut", "color"))
    expect_equal(res5, res6)

})
