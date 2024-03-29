test_that("replaceGreekSyms replaces Greek symbols", {
    expect_equal(replaceGreekSyms("β2-microglobulin", replace = "sym2word"),
                 "beta2-microglobulin")
    expect_equal(replaceGreekSyms("TCRγ/δ", replace = "sym2letter"), "TCRg/d")
    expect_equal(replaceGreekSyms("TCR gamma/delta", replace = "word2letter"),
                 "TCR g/d")
    expect_equal(replaceGreekSyms("β2-microglobulin"), "b2-microglobulin")
    expect_equal(replaceGreekSyms(c("IFN-γ R α chain", "FcεRIα")),
                                  c("IFN-g R a chain", "FceRIa"))
})
