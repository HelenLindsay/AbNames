# Tests for filter_by_union -----

data(diamonds, package = "ggplot2")
diamonds <- diamonds[1:20,] %>%
    dplyr::mutate(across(c("cut", "color", "clarity"), as.character))

exp_res1 <- diamonds %>% dplyr::filter(color == "I" | cut == "Fair")


test_that("filter_by_union selects the correct rows", {
    res1 <- filter_by_union(diamonds, data.frame(cut = "Fair", color = "I"))
    expect_equal(res1, exp_res1)

    # Value 4.05 appears in columns x and y.
    # filter_by_union should respect column
    res2 <- filter_by_union(diamonds, data.frame(x = 4.05))
    exp_res2 <- diamonds %>% dplyr::filter(x == 4.05)
    expect_equal(res2, exp_res2)

    # If one data.frame df and a row range is supplied, filter_by_union should
    # subset from df
    res3 <- filter_by_union(diamonds, rows = 10)
    exp_res3 <- diamonds %>%
        dplyr::filter(carat == 0.23 | cut == "Very Good" | color == "H" |
                          clarity == "VS1" | table == 61 )
    expect_equal(res3, exp_res3)

    # If two data.frames and a row.range are supplied, filter_by_union should
    # subset from second data.frame
    res4 <- filter_by_union(diamonds, diamonds[, c("cut", "color")], rows = 9:10)
    res5 <- filter_by_union(diamonds, diamonds[9:10, c("cut", "color")])
    expect_equal(res4, res5)

    # Using "by" to subset should be equivalent to explicitly subsetting df2
    res6 <- filter_by_union(diamonds, diamonds[9:10,], by = c("cut", "color"))
    expect_equal(res5, res6)
})


test_that("filter_by_union correctly handles names in 'by' argument", {
    res1 <- filter_by_union(diamonds, data.frame(Cut = "Fair", color = "I"),
                       by = c(cut = "Cut", "color"))
    expect_equal(res1, exp_res1)

    # Error if names of 'by' not in df
    expect_error(filter_by_union(diamonds, data.frame(Cut = "Fair", color = "I"),
                            by = c(fish = "Cut", "color")))

    # Error if values of 'by' not in df
    expect_error(filter_by_union(diamonds, data.frame(Cut = "Fair", color = "I"),
                            by = c(cut = "Cut", "fish")))
})


# Tests for group_by_any ----

test_that("group_by_any works", {
    df <- data.frame(A = c("a", "b", "c", "c", "d", "e"),
                     B = c("f", "g", "g", "h", "i", "j"),
                     C = c("k", "k", "l", "m", "n", "o"))
    group_by_any(df, c("A","B","C"))
})

# Tests for left_join_any ----

test_that("left_join_any updates correctly", {
    # df_y has A = a, D = k and A = a, D = j - both should match by col A
    # B = e has one match
    # A = c, B = f does not match, D should be NA
    # A = a, B = f - NA value from df_y should not be copied

    df_x <- data.frame(ID = c(1, 2, 3, 4),
                       A = c("a", "b", "c", "a"),
                       B = c("d", "e", "f", "f"),
                       C = c(NA, "h", "i", "f"))

    df_y <- data.frame(A = c("a", "a", "b", "a"),
                       B = c("j", "e", "b", NA),
                       C = c("g", "g", NA, NA),
                       D = c("k", "k", "k", "j"))

    exp <- data.frame(ID = c(1, 1, 2, 3, 4, 4),
                      A = c("a", "a", "b", "c", "a", "a"),
                      B = c("d", "d", "e", "f", "f", "f"),
                      C = c("g", "g", "h", "i", "f", "f"),
                      D = c("k", "j", "k", NA, "k", "j"))

    res <- left_join_any(df_x, df_y, cols = c("A","B")) %>% arrange(ID)
    expect_equal(res, exp)

    # To do:
    # Test without patching (no non-joining columns in common)
    # Test with patching
    # Test with updating

})
