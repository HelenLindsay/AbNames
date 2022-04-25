test_that(".tempColName returns a column name that doesn't already exist", {
    expect_equal(.tempColName(data.frame(x = 1)), "TEMP")
    expect_equal(.tempColName(data.frame(TEMP = 1)), "TEMP1")
    expect_equal(.tempColName(data.frame(x = 1), n = 3),
                 c("TEMP", "TEMP1", "TEMP2"))
    expect_equal(.tempColName(data.frame(TEMP = 1, TEMP2 = 2), n = 2),
                 c("TEMP1", "TEMP3"))
})
