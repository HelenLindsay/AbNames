df <- data.frame(Antigen = c("CD235a/b", "CD235ab", "CD66a/c/e",
                             "CD66ace", "HLA-A,B,C", "HLA-A/B/C",
                             "HLA-A2", "HLA-ABC", "HLA-DR",
                             "TCR alpha/beta", "TCR g/d", "TCRab"))

test_that(".separateSubunits works correctly with default patterns", {
    p1_res <- data.frame(Antigen = c("CD235a/b", "CD235a/b", "CD235ab",
            "CD235ab", "CD66a/c/e", "CD66a/c/e", "CD66a/c/e", "CD66ace",
            "CD66ace", "CD66ace", "HLA-A,B,C", "HLA-A/B/C", "HLA-A2", "HLA-ABC",
            "HLA-DR", "TCR alpha/beta", "TCR g/d", "TCR g/d", "TCRab", "TCRab"),
            subunit = c("CD235a", "CD235b", "CD235a", "CD235b", "CD66a",
            "CD66c", "CD66e", "CD66a", "CD66c", "CD66e", NA, NA, NA, NA,
            NA, NA, "TCRg", "TCRd", "TCRa", "TCRb"))


    #    "IL-3R", "CD235ab", "TCR alpha/beta",
    #    "TCRab", "HLA-DR", "B7-H2", "HLA-A,B,C", "HLA-A/B/C", "HLA-ABC",
    #    "DC-SIGN", "TCR gamma/delta", "TCRgd", "TCR g/d", "IL-2Rbeta",
    #    "DEC-205", "CD66a/c/e", "CSF-1R", "VE-cadherin", "IL-7Ra",
    #    "IL-7Ralpha", "IL-6Ra", "IL-4Ralpha", "IL-2Rb", "CD235a/b", "IL-21R",
    #    "CD66ace", "HLA-A2", "IL-10", "TCR a/b"))

    df <- data.frame(Antigen = c(""))
})
