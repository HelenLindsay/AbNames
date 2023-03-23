test_that("getCommonName input checks work", {
    expect_error(getCommonName(data.frame(Antigen_std = 1)))
    expect_error(getCommonName(data.frame(Antibody = 1)))
    expect_error(getCommonName(data.frame(Antigen = 1, n_matched = 2)))
    expect_error(getCommonName(data.frame(Antigen = 1), cols = "banana"))
})


test_that("getCommonName verbose works", {
    out <- capture_output(getCommonName(data.frame(Antigen = "CD4",
                                                   Clone = "CLONE_NM")))
    expect_true(grepl(out, "Grouping data using columns"))
})
