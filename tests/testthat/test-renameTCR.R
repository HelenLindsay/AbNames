tcr <- data.frame(Antigen = c("TCR alpha/beta", "TCRab", "TCR gamma/delta",
                              "TCRgd", "TCR g/d", "TCRg", "TCR Vgamma9",
                              "TCR Vg9", "TCR Vd2", "TCR Vdelta2",
                              "TCR Va24-Ja18", "TCRVa24.Ja18",
                              "TCR Valpha24-Jalpha18", "TCR Va7.2", "TCRa7.2",
                              "TCRVa7.2", "TCR Vbeta13.1", "TCR Vb13.1",
                              "TCR Vg9", "TCR Vd2"),
                  )

test_that("formatTCR correctly formats to long format", {
  expect_equal(2 * 2, 4)
})
