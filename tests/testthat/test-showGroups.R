test_that("showGroups correctly displays the chosen group", {
    df <- data.frame(x = rep(letters[1:3], each = 3))

    # Function should warn because df is ungrouped
    expect_warning(showGroups(df, interactive = FALSE))

    df <- dplyr::group_by(df, x)
    out <- capture_output(showGroups(df, i = 2, interactive = FALSE))
    expect_true(grepl("Group 2 of 3", out))

    # Expect and error as there are only 3 groups
    expect_error(showGroups(df, i = 4, interactive = FALSE))
})


test_that(".break_into_n returns the correct indices", {
    #.break_into_n returns the indices used by interactive function print_n

    # If piece size is > nrows, expect 1 piece
    expect_equal(.break_into_n(30, 40), list(starts = 1, ends = 30))
    # If piece size is a multiple of nrows, make sure end isn't missed
    expect_equal(.break_into_n(30, 10), list(starts = c(1,11,21),
                                             ends = c(10,20,30)))
    # Make sure last row isn't missed
    expect_equal(.break_into_n(30, 29),list(starts = c(1,30), ends = c(29,30)))
})
