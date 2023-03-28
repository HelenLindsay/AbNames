# .separateSubunits ----
# tests of individual splitting patterns
test_that(".separateSubunits works correctly with default patterns", {
    df <- data.frame(Antigen = c("CD235a/b", "CD235ab", "CD66a/c/e",
                                 "CD66ace", "HLA-A,B,C", "HLA-A/B/C",
                                 "HLA-A2", "HLA-ABC", "HLA-DR",
                                 "TCR alpha/beta", "TCR g/d", "TCRab"))

    p1 <- "^[A-Z0-9]{2,}[-\\. ]?([a-z\\/\\.]{2,6})$" # Only lowercase
    p2 <- "^[A-Z0-9]{2,}[-\\.]([A-Z\\/,\\.]{2,6})$" # Only uppercase

    p1_res <- c(rep(c("CD235a", "CD235b"), 2),
                rep(c("CD66a", "CD66c", "CD66e"), 2), rep(NA, 6),
                "TCRg", "TCRd", "TCRa", "TCRb")
    p2_res <- c(rep(NA, 4), rep(c("HLA-A", "HLA-B", "HLA-C"), 2),
                NA, sprintf("HLA-%s", c("A","B","C","D","R")), rep(NA, 3))

    df_p1_res <- .separateSubunits(df, "Antigen", "TMP3", p1,
                                   "%s%s", "TMP1", "TMP2")
    expect_equal(df_p1_res$TMP3, p1_res)

    df_p2_res <- .separateSubunits(df, "Antigen", "TMP4", p2, "%s-%s",
                                   "TMP1", "TMP2")
    expect_equal(df_p2_res$TMP4, p2_res)
})


# separateSubunits ----
test_that("separateSubunits doesn't split non-subunit genes", {
    not_subunits <- data.frame(
        Antigen = c("IL-3R", "B7-H2", "TCR gamma/delta",
                    "IL-2Rbeta", "DEC-205", "CSF-1R", "VE-cadherin",
                    "IL-7Ra", "IL-7Ralpha", "IL-6Ra", "IL-4Ralpha",
                    "IL-2Rb", "IL-21R", "HLA-A2", "IL-10"))
    res <- separateSubunits(not_subunits)
    # If none of the Antigens in "not_subunits" are split,
    # nrows shouldn't be changed
    expect_equal(nrow(res), nrow(not_subunits))
})
