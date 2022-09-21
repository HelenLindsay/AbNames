test_that("formatTCRv correctly formats to long format", {
    tcr <- c("TCR alpha/beta", "TCRab", "TCR g/d", "TCRg",
             "TCR Vgamma9", "TCR Vdelta2", "TCR Va24-Ja18", "TCRVa24.Ja18")
    n_sunits <- c(2,2,2,1,1,1,2,2)

    tcr_long <- c("T cell receptor alpha locus",
                  "T cell receptor beta locus",
                  "T cell receptor alpha locus",
                  "T cell receptor beta locus",
                  "T cell receptor gamma locus",
                  "T cell receptor delta locus",
                  "T cell receptor gamma locus",
                  "T cell receptor gamma variable 9",
                  "T cell receptor delta variable 2",
                  "T cell receptor alpha variable 24",
                  "T cell receptor alpha joining 18",
                  "T cell receptor alpha variable 24",
                  "T cell receptor alpha joining 18")

    expect_equal(as.data.frame(formatTCRv(tcr)),
                 data.frame(TCR = rep(tcr, n_sunits), TCR_long = tcr_long))
})


test_that("formatTCR correctly merges results", {
    tcr <- c("TCR alpha/beta", "TCRab", "TCR g/d", "TCRg",
             "TCR Vgamma9", "TCR Vdelta2", "TCR Va24-Ja18", "TCRVa24.Ja18")

    df <- data.frame(ID = c("A", "B"), tcr = tcr[1:2])
    res <- data.frame(ID = rep(c("A", "B"), each = 2),
                      tcr = rep(tcr[1:2], each = 2),
                      long_name = rep(sprintf("T cell receptor %s locus",
                                              c("alpha", "beta")), 2))
    expect_equal(formatTCR(df, "tcr", "long_name"), res)

})
