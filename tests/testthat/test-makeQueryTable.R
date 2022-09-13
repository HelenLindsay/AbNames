# Collecting some cases for testing
#test <- data.frame(ID = 1:7,
#               Antigen = c("CD274.1", "CD274_PROT",
#                           "CD274_B7-H1_PD-L1",
#                           "anti-human_CD274_(B7-H1_and_PD-L1)",
#                           "CD274|CD274|AHS0004|pAbO",
#                                     "CD274 (B7-H1, PD-L1)",
#                           "PDL1 (CD274)"))


test_that("makeQueryTable correctly handles brackets and commas", {
    query_df <- data.frame(ID = "test",
                           Antigen = "CD274 (B7-H1, PD-L1)")

    qdf <- makeQueryTable(query_df)
    should_include <- c("CD274", "B7-H1", "PD-L1")
    expect_equal(intersect(qdf$value, should_include), should_include)
})


test_that("upperSquish correctly ignores groups of numbers", {
    test <- c("CD3.1", "CD3-A", "Nectin-2", "IL-2Rb", "CD8a", "LFA.3")
    exp <- c(NA, "CD3A", "NECTIN2", "IL2RB", "CD8A", "LFA3")
    expect_equal(upperSquish(test), exp)
})
