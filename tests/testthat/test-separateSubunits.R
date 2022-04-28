test_that("separateSubunits", {

    df <- data.frame(Antigen = c("CD235a/b", "CD235ab", "CD66a/c/e",
                                 "CD66ace", "HLA-A,B,C", "HLA-A/B/C",
                                 "HLA-A2", "HLA-ABC", "HLA-DR",
                                 "TCR alpha/beta", "TCR g/d", "TCRab"))

    #    "IL-3R", "CD235ab", "TCR alpha/beta",
    #    "TCRab", "HLA-DR", "B7-H2", "HLA-A,B,C", "HLA-A/B/C", "HLA-ABC",
    #    "DC-SIGN", "TCR gamma/delta", "TCRgd", "TCR g/d", "IL-2Rbeta",
    #    "DEC-205", "CD66a/c/e", "CSF-1R", "VE-cadherin", "IL-7Ra",
    #    "IL-7Ralpha", "IL-6Ra", "IL-4Ralpha", "IL-2Rb", "CD235a/b", "IL-21R",
    #    "CD66ace", "HLA-A2", "IL-10", "TCR a/b"))

    df <- data.frame(Antigen = c(""))
})
