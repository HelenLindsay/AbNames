df <- data.frame(A = rep(1:3, each = 2),
                 B = rep(4:6, 2),
                 C = c(7,7,10,8,10,10),
                 nB = c(11:16))

test_that(".addNPerGroup works", {
    exp <- cbind(df, data.frame(nB = 2,
                                 nC = rep(c(1,2,1), each = 2)))

    res1 <- .addNPerGroup(df, A, c("B", "C"))

    expect_equal(as.data.frame(res1), exp1)
})
