
test_that("searchAliases does not accept partial subunit matches", {
    qdf <- data.frame(ID = rep(c("CD66acz", "CD66ace"), each = 3),
                      value = c("CD66a", "CD66c", "CD66z",
                                "CD66a", "CD66c", "CD66e"),
                      name = rep("subunit", 6))

    res <- searchAliases(qdf)
    # Expect that a match is found only for "CD66ace", not for made up "CD66acz"
    expect_equal(unique(res$ID), "CD66ace")
})


test_that("searchAliases prefers manual symbols to HGNC ID or symbol", {

    # These will both HGNC ID HGNC:9666
    # CD45R (manual lookup) should be preferred over PTP
    qdf <- data.frame(ID = rep("CD45R_PTPRC", 2),
                      value = c("CD45R", "PTPRC"),
                      name = rep(NA, 2))
    res <- searchAliases(qdf)
    expect_equal(res$SOURCE, "MANUAL_LOOKUP")
})


test_that("searchAliases prefers hgnc symbols to previous symbols", {
    # HGNC_SYMBOL is CD8A, CD8 is a PREVIOUS_SYMBOL
    qdf <- data.frame(ID = rep("CD8_CD8A", 2),
                      value = c("CD8", "CD8A"),
                      name = NA)
    res <- searchAliases(qdf)
    expect_equal(res$symbol_type, "HGNC_SYMBOL")
})


test_that("searchAliases aggregates symbols if no official symbol is found", {
    # CD227 and BT3.1 map to different genes (deliberate typo for CD277)
    # CD134 and OX40 are both aliases but neither is the official symbol

    # Here we expect the first antigen to map to 2 different genes,
    # and the second antigen to map to the same gene via 2 different symbols

    qdf <- data.frame(ID = rep(c("CD277_BT3.1", "CD134_OX40"), each = 2),
                      value = c("CD227", "BT3.1",
                                "CD134", "OX40"),
                      name = NA)

    res <- searchAliases(qdf)
    expect_equal(res$value, c("CD227", "BT3.1", "CD134|OX40"))
})
