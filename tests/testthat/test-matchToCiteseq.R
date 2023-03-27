test_that("matchToCiteseq input and output checks work", {
    # Input data.frame must have a column "Antigen"
    expect_error(matchToCiteseq(data.frame(Antibody = "CD14")))
    # Warn if "cols" are not present in citeseq
    expect_warning(matchToCiteseq(data.frame(Antigen = "CD14"),
                                  cols = "Antibody"))
})


test_that("matchToCiteseq correctly calls group_by_any", {

    custom_abs <- data.frame(ID = LETTERS[1:8],
                             Antigen = c("CD27", "CD28", "CD38", "CD56",
                                         "CD57", "ICOS", "CD27", "CD38(EXTRA)"),
                             Cat_Number = rep("custom made", 8),
                             Clone = c("REA499", "CD28.2", "REA572", "REA196",
                                       "QA17A04", "C398.4A", "NONSENSE",
                                       "REA572"))
    res <- matchToCiteseq(custom_abs)
    # These entries shouldn't be matched by Cat_Number = "custom made",
    # so we expect that col "Antigen_std" should have 6 unique values
    expect_equal(length(unique(res$Antigen_std)), 6)
})
