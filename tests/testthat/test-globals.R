test_that("quick tests that global functions run without error",{
    expect_true(class(GREEKSYM_TO_NAME()) == "character")
    expect_true(class(GREEKSYM_TO_LETTER()) == "character")
    expect_true(class(GREEKNAME_TO_LETTER()) == "character")
    expect_true(class(CLONE_DUPS()) == "data.frame")
    expect_true(class(ISOTYPE_CONTROLS()) == "data.frame")
})
