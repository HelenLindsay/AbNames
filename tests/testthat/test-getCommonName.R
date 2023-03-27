test_that("getCommonName input checks work", {
    expect_error(getCommonName(data.frame(Antigen_std = 1)))
    expect_error(getCommonName(data.frame(Antibody = 1)))
    expect_error(getCommonName(data.frame(Antigen = 1, n_matched = 2)))
    expect_error(getCommonName(data.frame(Antigen = 1), cols = "banana"))
})

test_that("getCommonName verbose works", {
    expect_message(getCommonName(data.frame(Antigen = "CD4",
                                            Clone = "CLONE_NM")),
                   "Grouping data using columns")
})

test_that("getCommonName and group_by_any correctly handle ignore arg", {
    # These match by clone but should be grouped into three groups because
    # only the clone is shared for the CD197s
    df <- data.frame(ID = c("CD197__Mair_2020", "CD197 (CCR7)__Trzupek_2021",
                            "CD197 (CCR7)__Trzupek_2020", "HLA-E__LeCoz_2021",
                            "HLA-E__Frangieh_2021"),
                     Antigen = c("CD197", "CD197 (CCR7)", "CD197 (CCR7)",
                                 "HLA-E", "HLA-E"),
                     Cat_Number = c(NA, "940014", NA, "342619", NA),
                     Oligo_ID = c("AHS0007", NA, NA, "0918", "0918"),
                     Clone = c("3D12", "3D12", "3D12", "3D12", "3D12"))

    res <- getCommonName(df)
    expect_equal(length(unique(res$group)), 3)

    # custom antibodies including one pair that matches by Antigen
    # and one pair that matches by Clone
    custom_abs <- data.frame(Antigen = c("CD27", "CD28", "CD38", "CD56",
                                         "CD57", "ICOS", "CD27", "CD38(EXTRA)"),
                             Cat_Number = rep("custom made", 8),
                             Clone = c("REA499", "CD28.2", "REA572", "REA196",
                                       "QA17A04", "C398.4A", "NONSENSE",
                                       "REA572"))
    res <- getCommonName(custom_abs)
    expect_equal(length(unique(res$group)), 6)
})


test_that("getCommonName does not group by NA", {
    df <- data.frame(Antigen = c(rep(c("B220 (CD45R)", "B7-H4"), each = 2),
                                 "C-kit (CD117)", "C-KIT"),
                     Cat_Number = rep(NA_character_, 6),
                     HGNC_ID = rep(c("HGNC:9666", "HGNC:28873", "HGNC:6342"),
                                   each = 2))

    # Without using HGNC_ID, there should be 4 groups
    res <- getCommonName(df, cols = c("Antigen", "Cat_Number"))
    expect_equal(length(unique(res$group)), 4)

    # If HGNC_ID is also used, there should now be 3 groups
    res <- getCommonName(df, cols = colnames(df))
    expect_equal(length(unique(res$group)), 3)
})
