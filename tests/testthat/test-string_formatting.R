
# addID ----

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

# splitUnnest ----

test_that("splitUnnest correctly splits at brackets", {
    df <- data.frame(Antigen = c("CD279 (PD-1)", "CD3", "Mac-2 (Galectin-3)"))
    exp_res <- data.frame(Antigen = c(rep("CD279 (PD-1)", 2),"CD3",
                                      rep("Mac-2 (Galectin-3)", 2)),
                          Split = c("CD279", "PD-1", "CD3",
                                    "Mac-2", "Galectin-3"))

    # Test with adding a new column
    expect_equal(data.frame(splitUnnest(df, "Antigen", new_col = "Split")),
                 exp_res)

    # Test overwriting the existing column
    expect_equal(data.frame(splitUnnest(df, "Antigen")),
                 data.frame(Antigen = exp_res$Split))
})


test_that("splitUnnest works with exclude", {
    ag <- c("CD8a/CD8A", "TCR alpha/beta", "erbB2/HER-2")
    exp_res <- data.frame(Antigen = c(rep(ag[1], 2), ag[2], rep(ag[3], 2)),
                      Split = c("CD8a", "CD8A", "TCR alpha/beta",
                                "erbB2", "HER-2"))
    res <- splitUnnest(data.frame(Antigen = ag), split = "\\/",
                       exclude = "TCR", new_col = "Split")
    expect_equal(data.frame(res), exp_res)
})


# gsubAb ----



