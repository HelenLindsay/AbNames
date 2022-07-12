# CD158f/ KIR2DL5
# https://doi.org/10.3389/fimmu.2012.00289
# KIR2DL5 is highly polymorphic and exhibits copy number variation
# KIR2DL5A (telomeric) and KIR2DL5B (centromeric, usually inactivated)
# KIR2DL5 is represented in the Immuno Polymorphism Database by 15 KIR2DL5A
# and 25 KIR2DL5B allels.  One single polymorphism distinguishes A and B

IFN_gamma <- paste("https://proconsortium.org/app/entry/PR:000001361/",
                   "IFN-gamma receptor 1 and IFN-gamma-R-alpha are synonyms")

HLA_DR <- paste("Grep in protein ontology for HLA-DR gives",
                "wrong results PR:P04233, PR:000001822")

CD66 <- paste0("https://www.bdbiosciences.com/en-ch/products/reagents/",
               "flow-cytometry-reagents/research-reagents/",
               "single-color-antibodies-ruo/pe-mouse-anti-human-cd66.551480")


manual_matches <- tibble::tribble(
    ~Antigen, ~Clone, ~HGNC_Symbol, ~HGNC_ID, ~PRO_ID, ~Comments, ~Cat_Number,
    "Annexin V", "", "ANXA5", "HGNC:543", "", "", "",
    "Mac-2", "", "LGALS3", "HGNC:6563", "", "", "",
    "FAS.L", "", "FASLG", "HGNC:11936", "", "", "",
    "IFN-g R a chain", "", "IFNGR1", "HGNC:5439", "PR:000001361", "IFN_gamma", "",
    "TCR Va7.2", "","TRAV7", "HGNC:12145", "", "", "",
    "cKIT", "","KIT", "HGNC:6342", "", "", "",
    "CD77", "","A4GALT", "HGNC:18149", "",
        "https://www.sinobiological.com/research/cd-antigens/cd77", "",
    "Podocalyxin", "","PODXL", "HGNC:9171", "",
        "Podocalyxin is podocalxyin-like", "",
    "[Ff]olate [Rr]eceptor b", "FOLR2", "HGNC:3793", "", "", "", "",
    "PE", "", NA, NA, NA,
        "anti-phycoerythrin, for binding PE-antibody labeled cells", "",
    "integrin b7", "","ITGB7", "HGNC:6162", "", "","",
    "HLA-DR", "", "", "","PR:P01903, PR:000050459", "", "",
    "CD11a/CD18", "", "ITGB2, ITGAL", "", "", "", "",
    "CD66","B1.1", "CEACAM1, CEACAM6, CEACAM3, CEACAM5", "", "", "", "",
    "Tau", "", "MAPT", "HGNC:6893", "PR:000027448, PR:000027447",
        "phospho-tau (thr181)", ""
)
