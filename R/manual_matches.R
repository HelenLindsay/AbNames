IFN_gamma <- paste("https://proconsortium.org/app/entry/PR:000001361/",
                   "IFN-gamma receptor 1 and IFN-gamma-R-alpha are synonyms")

manual_matches <- tibble::tribble(
    ~Antigen, ~HGNC_Symbol, ~HGNC_ID, ~PRO_ID, ~Comments, ~Cat_Number,
    "Annexin V", "ANXA5", "HGNC:543", "", "", "",
    "Mac-2", "LGALS3", "HGNC:6563", "", "", "",
    "FAS.L", "FASLG", "HGNC:11936", "", "", "",
    "IFN-g R a chain", "IFNGR1", "HGNC:5439", "PR:000001361", IFN_gamma, "",
    "TCR Va7.2", "TRAV7", "HGNC:12145", "", "", "",
    "cKIT", "KIT", "HGNC:6342", "", "", "",
    "CD77", "A4GALT", "HGNC:18149", "",
        "https://www.sinobiological.com/research/cd-antigens/cd77", "",
    "Podocalyxin", "PODXL", "HGNC:9171", "",
        "Podocalyxin is podocalxyin-like", "",
    "[Ff]olate [Rr]eceptor b", "FOLR2", "HGNC:3793", "", "", "",
    "PE", NA, NA, NA,
        "anti-phycoerythrin, for binding PE-antibody labeled cells", "",
    "integrin b7", "ITGB7", "HGNC:6162", "", "",""
)
