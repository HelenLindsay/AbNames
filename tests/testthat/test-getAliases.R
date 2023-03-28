test_that("getAliases correctly switches between ALT_ID and HGNC_ID", {

    # There should be no aliases of CD45RA using ALT_ID, as this was
    # manually added
    expect_equal(getAliases("CD45RA")$value, "CD45RA")

    # If using HGNC_ID, we expect that CD45 and CD45RO and B220 are aliases
    exp_aliases <- c("CD45RA", "CD45RO", "CD45", "B220")
    expect_true(all(exp_aliases %in%
                        getAliases("CD45RA", by = "HGNC_ID")$value))
})

test_that("getAliases prints a message if antigen is not found", {
    expect_message(getAliases("BANANA"), "not found")
    expect_no_message(getAliases("BANANA", verbose = FALSE))
})
