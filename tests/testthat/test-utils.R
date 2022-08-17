# Tests for .tempColName -----

test_that(".tempColName returns a column name that doesn't already exist", {
    expect_equal(.tempColName(data.frame(x = 1)), "TEMP")
    expect_equal(.tempColName(data.frame(TEMP = 1)), "TEMP1")
    expect_equal(.tempColName(data.frame(x = 1), n = 3),
                 c("TEMP", "TEMP1", "TEMP2"))
    expect_equal(.tempColName(data.frame(TEMP = 1, TEMP2 = 2), n = 2),
                 c("TEMP1", "TEMP3"))
})


# Tests for .wordGrep -----

test_that(".wordGrep handles pattern at all positions", {
    expect_equal(.wordGrep("CD3", "CD33"), FALSE)
    expect_equal(.wordGrep("CD3", "in the CD3 middle"), TRUE)
    expect_equal(.wordGrep("CD3", "sentence about CD3."), TRUE)
    expect_equal(.wordGrep("CD3", "two CD3. Sentences."), TRUE)
    expect_equal(.wordGrep("CD3", "CD3"), TRUE)
    expect_equal(.wordGrep("CD3", "CD3.A"), FALSE)
    expect_equal(.wordGrep("CD3", "Sentence with CD3.A"), FALSE)
})


# Tests for .stopIfColExists -----

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
