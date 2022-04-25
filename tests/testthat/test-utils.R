test_that(".tempColName returns a column name that doesn't already exist", {
    expect_equal(.tempColName(data.frame(x = 1)), "TEMP")
    expect_equal(.tempColName(data.frame(TEMP = 1)), "TEMP1")
})
