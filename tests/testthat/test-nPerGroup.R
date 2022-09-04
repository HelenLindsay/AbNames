df <- data.frame(A = rep(1:3, each = 2),
                 B = rep(4:6, 2),
                 C = c(7,7,10,8,10,10),
                 D = rep(11:12, c(4,2)),
                 nB = c(11:16))

test_that(".addNPerGroup works with multiple columns", {
    # One grouping column, two columns for n_distinct
    exp1 <- cbind(df[, c("A","B","C","D")], data.frame(nB = 2,
                                 nC = rep(c(1,2,1), each = 2)))

    # Both quoted
    res1_3 <- .addNPerGroup(df, "A", c("B", "C"))
    expect_equal(as.data.frame(res1_3), exp1)


    # Two grouping columns, one column for n_distinct
    exp2 <- cbind(df[,c("A","B","C","D")],
                 data.frame(nB = c(2,2,1,1,2,2)))
    res2_1 <- .addNPerGroup(df, c("A","C"), "B")
    expect_equal(as.data.frame(res2_1), exp2)

    # Two grouping columns, two columns for n_distinct
    res3 <- .addNPerGroup(df, c("A","C"), c("B", "D"))
    expect_equal(as.data.frame(res3), cbind(exp2, data.frame(nD = 1)))
})


test_that("nPerGroup works correctly with incomplete cases", {
    df[2, 1] <- NA
    df[3, 2] <- NA
    expect_error(nPerGroup(df, c("A","C"), c("B", "D")),
                 regex = " nB |Please")
    df$nB <- NULL

    res <- nPerGroup(df, c("A","C"), c("B", "D"))
    # NA columns are moved to the bottom
    exp_res <- cbind(df[c(1,3:6,2),], data.frame(nB = c(1, 0, 1, 2, 2, NA),
                                                 nD = c(1, 1, 1, 1, 1, NA)))
    rownames(exp_res) <- NULL
    expect_equal(as.data.frame(res), exp_res)
})
