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
