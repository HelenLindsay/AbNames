test_that("searchTotalseq searches by Cat_Number or Antigen", {
    cat_num <- c("306623", "124235", "331609")
    antigen <- c("CD235ab", "CD27", "CD276")

    # Given only the Cat_Number, the matching antigen should be found
    res_cat_number = searchTotalseq(data.frame(Cat_Number = cat_num))
    expect_equal(res_cat_number$Antigen, antigen)

    # Given the Antigen, the ENSEMBL_ID should still match the results above
    expect_equal(searchTotalseq(data.frame(Antigen = antigen))$ENSEMBL_ID,
                 res_cat_number$ENSEMBL_ID)
})


test_that("searchTotalseq matches by multiple columns and fills Antigen", {
    antigen = c("CD235ab", "CD27", "CD276")
    dat <- data.frame(Cat_Number = c("306623", NA, "331609"),
                      Antigen = c(NA, "CD27", "CD276"))

    # Note that left_join_any does sequential matching so row order can change
    expect_setequal(searchTotalseq(dat)$Antigen, antigen)
})


test_that("searchTotalseq doesn't overwrite existing identifier", {
    dat <- data.frame(Cat_Number = c("306623", NA, "331609"),
                      Antigen = c(NA, "CD27", NA),
                      ENSEMBL_ID = c(NA, NA, "ENSG00000103855"))

    expect_equal(searchTotalseq(dat)$Antigen, c("CD235ab", "CD27", NA))
})


test_that("searchTotalseq doesn't fill Cat_Number", {
    antigen <- c("CD235ab", "CD27", "CD276")
    # NA Cat_Number should not be filled
    df <- data.frame(Antigen = antigen, Cat_Number = NA_character_)
    expect_equal(searchTotalseq(df)$Cat_Number, df$Cat_Number)
})
