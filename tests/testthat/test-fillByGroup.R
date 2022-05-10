# Data for use in tests ----

df <- data.frame(A = c(rep("A", 3), NA, rep("B", 4), rep("C", 3), NA),
                 B = c(rep("A", 2), NA, "A", rep("B", 4), rep("C", 3), NA),
                 C = c(NA, 1, 1, 1, 2, 3, NA, 2, 4, 4, NA, 5))

df_c <- df %>% dplyr::filter(! is.na(A) & ! is.na(B))
df_d <- c(6, NA, 6, 7, 8, NA, 8, NA, 9, 8, NA, 10)

# Simple filling function using tidyr::fill
partial_f <- function(df, col, gp){
    df %>%
        dplyr::group_by(!!!syms(gp)) %>%
        tidyr::fill(!!sym(col), .direction = "updown") %>%
        dplyr::ungroup()
}

test_that("groupMode correctly adds a new column if required", {
    res_col <- c(1,1,2,3,2,2,4,4,4)
    expect_equal(as.data.frame(groupMode(df_c, "C", c("A","B"))),
                 data.frame(A = df_c$A, B = df_c$B, C = res_col))
    expect_equal(as.data.frame(groupMode(df_c, "C", c("A","B"), "D")),
                 cbind(df_c, D = res_col))
})


test_that("groupMode checks minimum count per group", {
    res_col <- c(1,1,2,3,NA,2,4,4,4)
    expect_equal(as.data.frame(groupMode(df_c, "C", c("A","B"), min_n = 3)),
                 df_c)
})


test_that("fillByGroup with option=stop stops with multiple values", {
    expect_error(fillByGroup(df, c("A", "B"), "C", multiple = "stop"))
})


test_that("fillByGroup with option=mode fills multiple values", {
    exp_res <- data.frame(A = c(rep("A", 2), rep("B", 4),
                                rep("C", 3), "A", NA, NA),
                          B = c(rep("A", 2), rep("B", 4),
                                rep("C", 3), NA, "A", NA),
                     C = c(1, 1, 2, 3, 2, 2, 4, 4, 4, 1, 1, 5))
    res <- fillByGroup(df, group = c("A", "B"), fill = "C", multiple = "mode")
    expect_equal(as.data.frame(res), exp_res)
})


test_that(".freducePartial correctly handles ellipsis", {
    # Applying .freducePartial (in a loop) should equal applying
    # fill directly (vectorised)

    df$D <- df_d
    res <- .freducePartial(df, partial_f, cls = "col",
                           gp = c("A", "B"), col = c("C", "D"))
    exp_res <- df %>%
        dplyr::group_by(A, B) %>%
        tidyr::fill(C, D, .direction = "updown")

    expect_equal(as.data.frame(res), as.data.frame(exp_res))
})


test_that("fillByGroup can fill majority value for multiple columns", {
    df$D <- df_d
    res <- fillByGroup(df, group = c("A", "B"),
                          fill = c("C", "D"), multiple = "mode")
    exp_res <- data.frame(
        A = c("A", "A", "B", "B", "B", "B", "C", "C", "C", "A", NA, NA),
        B = c("A", "A", "B", "B", "B", "B", "C", "C", "C", NA, "A", NA),
        C = c(1, 1, 2, 3, 2, 2, 4, 4, 4, 1, 1, 5),
        D = c(6, 6, 8, 8, 8, 8, 9, 8, 9, 6, 7, 10))
    expect_equal(as.data.frame(res), exp_res)
})

