test_that("matchToCiteseq input and output checks work", {
    # Input data.frame must have a column "Antigen"
    expect_error(matchToCiteseq(data.frame(Antibody = "CD14")))
    # Warn if "cols" are not present in citeseq
    expect_warning(matchToCiteseq(data.frame(Antigen = "CD14"),
                                  cols = "Antibody"))
    # Expect warning when incorrect Cat_Number used
    expect_warning(matchToCiteseq(data.frame(Antigen = "CD14",
                                             Cat_Number = "329527")))
})
