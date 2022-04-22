


test_that("existing columns are not overwritten", {
    expect_error(addID(data.frame(x = 1), id_cols = "dummy", new_col = "x"))
})

test_that("addID warns if ID is not unique", {
    df <- data.frame(x = rep("A", 2), y = rep("B", 2))
    expect_warning(addID(df, id_cols = c("x", "y")))
})

test_that("addID allows more than two ID columns", {
    df <- data.frame(x = rep("A", 2), y = rep("B", 2), z = c("C","D"))
    result <- addID(df, c("x","y","z"))$ID
    expect_equal(result, c("A-B-C", "A-B-D"))
})
