
test_that("groupMode correctly adds a new column if required", {
    dat <- fillByGroup_data()
    res_col <- c(1,1,2,3,2,2,4,4,4)
    expect_equal(as.data.frame(groupMode(dat$df_c, "C", c("A","B"))),
                 data.frame(A = dat$df_c$A, B = dat$df_c$B, C = res_col))

    # If a new column is provided, groupMode now always takes the mode
    overwrite_res_col <- c(1,1,2,2,2,2,4,4,4)
    expect_equal(as.data.frame(groupMode(dat$df_c, "C", c("A","B"), "D")),
                 cbind(dat$df_c, D = overwrite_res_col))
})


test_that("groupMode checks minimum count per group", {
    dat <- fillByGroup_data()
    res_col <- c(1,1,2,3,NA,2,4,4,4)
    expect_equal(as.data.frame(groupMode(dat$df_c, "C", c("A","B"), min_n = 3)),
                 dat$df_c)
})


test_that("groupMode can correctly add a new column", {
    df <- data.frame(Antigen=c("CD274", "CD274", "PD-L1"),
                     Clone=rep("29E.2A3", 3))
    res <- groupMode(df, cl="Antigen", gp="Clone", new_cl="Antigen_mode")
    expect_equal(res$Antigen_mode, rep("CD274", 3))
})


test_that("fillByGroup with option=stop stops with multiple values", {
    dat <- fillByGroup_data()
    expect_error(fillByGroup(dat$df, c("A", "B"), "C", multiple = "stop"))
})


test_that("fillByGroup with option=mode fills multiple values", {
    dat <- fillByGroup_data()
    exp_res <- data.frame(A = c(rep("A", 2), rep("B", 4),
                                rep("C", 3), "A", NA, NA),
                          B = c(rep("A", 2), rep("B", 4),
                                rep("C", 3), NA, "A", NA),
                     C = c(1, 1, 2, 3, 2, 2, 4, 4, 4, 1, 1, 5))
    res <- fillByGroup(dat$df, group = c("A", "B"),
                       fill = "C", multiple = "mode")
    expect_equal(as.data.frame(res), exp_res)
})


test_that(".freducePartial correctly handles ellipsis", {
    # Applying .freducePartial (in a loop) should equal applying
    # fill directly (vectorised)
    dat <- fillByGroup_data()

    dat$df$D <- dat$df_d
    res <- .freducePartial(dat$df, dat$partial_f, cls = "col",
                           gp = c("A", "B"), col = c("C", "D"))
    exp_res <- dat$df %>%
        dplyr::group_by(A, B) %>%
        tidyr::fill(C, D, .direction = "updown")

    expect_equal(as.data.frame(res), as.data.frame(exp_res))
})


test_that("fillByGroup can fill majority value for multiple columns", {
    dat <- fillByGroup_data()

    dat$df$D <- dat$df_d
    res <- fillByGroup(dat$df, group = c("A", "B"),
                       fill = c("C", "D"), multiple = "mode")
    exp_res <- data.frame(
        A = c("A", "A", "B", "B", "B", "B", "C", "C", "C", "A", NA, NA),
        B = c("A", "A", "B", "B", "B", "B", "C", "C", "C", NA, "A", NA),
        C = c(1, 1, 2, 3, 2, 2, 4, 4, 4, 1, 1, 5),
        D = c(6, 6, 8, 8, 8, 8, 9, 8, 9, 6, 7, 10))
    expect_equal(as.data.frame(res), exp_res)
})


test_that("fillByGroup can ignore multiple modes", {
    # In the group A = "C" and B = "C", column "D" has 2 possible values,
    # one entry each
    dat <- fillByGroup_data()

    dat$df$D <- dat$df_d
    res <- fillByGroup(dat$df, group = c("A", "B"),
                       fill = c("C", "D"), multiple = "ignore")

    exp_res <- data.frame(
        A = c("A", "A", "B", "B", "B", "B", "C", "C", "C", "A", NA, NA),
        B = c("A", "A", "B", "B", "B", "B", "C", "C", "C", NA, "A", NA),
        C = c(1, 1, 2, 3, 2, 2, 4, 4, 4, 1, 1, 5),
        D = c(6, 6, 8, 8, 8, 8, 9, 8, NA, 6, 7, 10))
    expect_equal(as.data.frame(res), exp_res)
})


test_that("fillByGroup correctly overwrites if specified", {
    # In the group A = "C" and B = "C", column "D" has 2 possible values,
    # one entry each - these should stay NA, col C with a majority value
    # should be overwritten
    dat <- fillByGroup_data()

    dat$df$D <- dat$df_d
    res <- fillByGroup(dat$df, group = c("A", "B"), fill = c("C", "D"),
                       multiple = "ignore", method = "all")

    exp_res <- data.frame(
        A = c("A", "A", "B", "B", "B", "B", "C", "C", "C", "A", NA, NA),
        B = c("A", "A", "B", "B", "B", "B", "C", "C", "C", NA, "A", NA),
        C = c(1, 1, 2, 2, 2, 2, 4, 4, 4, 1, 1, 5),
        D = c(6, 6, 8, 8, 8, 8, NA, NA, NA, 6, 7, 10))
    expect_equal(as.data.frame(res), exp_res)
})
