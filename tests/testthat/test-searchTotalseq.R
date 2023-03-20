test_that("searchTotalseq searches by Cat_Number or Antigen", {
    cat_num <- c("306623", "124235", "331609")
    antigen <- c("CD235ab", "CD27", "CD276")

    expect_equal(searchTotalseq(data.frame(Cat_Number = cat_num))$Antigen,
                 antigen)

    # Cat_Number only returned if requested
    df <- data.frame(Antigen = antigen, Cat_Number = NA_character_)
    expect_equal(searchTotalseq(df)$Cat_Number, cat_num)
})

test_that("searchTotalseq matches by multiple columns", {
    dat <- data.frame(Cat_Num = c("306623", NA),
                      Antigen = c())

    cat_num <- c("306623", "124235", "331609")
    antigen <- c("CD235ab", "CD27", "CD276")
    expect_equal(searchTotalseq(data.frame(Cat_Number = cat_num))$Antigen,
                 antigen)
    expect_equal(searchTotalseq(data.frame(Cat_Number = cat_num))$Cat_Number,
                 cat_num)
})


test_that("searchTotalseq doesn't overwrite existing identifier", {
    dat <- data.frame(Cat_Num = c("306623", NA, "331609"),
                      Antigen = c(NA, "CD27", "CD276"),
                      ENSEMBL_ID = c(NA, NA, "ENSG00000103855"))

    # But maybe it should warn?
})
