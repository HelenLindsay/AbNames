# CD158f/ KIR2DL5
# https://doi.org/10.3389/fimmu.2012.00289
# KIR2DL5 is highly polymorphic and exhibits copy number variation
# KIR2DL5A (telomeric) and KIR2DL5B (centromeric, usually inactivated)
# KIR2DL5 is represented in the Immuno Polymorphism Database by 15 KIR2DL5A
# and 25 KIR2DL5B alleles.  One single polymorphism distinguishes A and B

IFN_gamma <- paste("https://proconsortium.org/app/entry/PR:000001361/",
                   "IFN-gamma receptor 1 and IFN-gamma-R-alpha are synonyms")

HLA_DR <- paste("Grep in protein ontology for HLA-DR gives",
                "wrong results PR:P04233, PR:000001822")

HLA_ABC <- paste("Biolegend website: Clone W6/32 recognizes ",
                 "human beta2-microglobulin")

CD66 <- paste0("https://www.bdbiosciences.com/en-ch/products/reagents/",
               "flow-cytometry-reagents/research-reagents/",
               "single-color-antibodies-ruo/pe-mouse-anti-human-cd66.551480")

#CD158 <- paste0("From BioLegend website: ",
#                "mAb HP-MA4 reacts with KIR2DL1 (CD158a), ",
#                "KIR2DS1 (CD158h), KIR2DS3, and KIR2DS5 (CD158g).")
#IgG <- paste0("From BioLegend website: clone M1310G05 has stronger affinity ",
#              "for IgG1 and IgG3 than for IgG2 and IgG4")

TRAV7 <- paste0("Vα7.2 TCR with Jα33 forms invariant T cell receptor")

CLA <- paste0("From BioLegend website: ",
              "CLA is a scarbohydrate epitope of sialic acid ",
              "and fucose-modified P-selectin glycoprotein ligand-1 (PSGL-1) ",
              "Ligand for E-selectin, P-selectin, and L-selectin.")


# Use ALT_ID when it is a modification, a complex, or not a protein
manual_matches <- tibble::tribble(~Antigen, ~Clone, ~HGNC_Symbol, ~HGNC_ID,
                                  ~PRO_ID, ~ALT_ID, ~Comments,

    #"Tau Phospho (Thr181)", "M7004D06", "MAPT", "HGNC:6893", "PR:000027448,
    #    PR:000027447", "PR:000027448",
    #    "PRO-short-label: hMAPT/iso:Tau-F/Phos:1",

#    "Mac-2", "M3/38", "LGALS3", "HGNC:6563", "", "", "", "",


    # This can be matched via the clone
    "FAS.L", "NOK-1", "FASLG", "HGNC:11936", "", "", "", "",

    "IFN-g R a chain", "GIR-208", "IFNGR1", "HGNC:5439",
        "PR:000001361", "", IFN_gamma,

    "TCR Va7.2", "3C10","TRAV7", "HGNC:12145", "", TRAV7,
    #"TCR Vb13.1", "H131", "", "", "", "", "",

    "cKIT", "104D2","KIT", "HGNC:6342", "", "", "",

    "CD77", "","A4GALT", "HGNC:18149", "", "",
        "https://www.sinobiological.com/research/cd-antigens/cd77",

    "[Ff]olate [Rr]eceptor b", "FOLR2", "HGNC:3793", "", "", "", "",

   # "PE", "PE001", NA, NA, NA, "",
   #     "anti-phycoerythrin, for binding PE-antibody labeled cells",

    "integrin b7", "FIB504","ITGB7", "HGNC:6162", "", "", "",

    "HLA-DR", "L243", NA, NA, "PR:P01903, PR:000050459", "", HLA_DR,
    "HLA.A.B.C", "W6/32", "B2M", "HGNC:914", "", "", HLA_ABC,

    "CD11a/CD18", "", "ITGB2, ITGAL", "", "", "", "",

   # "CD66","B1.1", "CEACAM1, CEACAM6, CEACAM3, CEACAM5", "", "", "", CD66,




    "CD158", "HP-MA4", "KIR2DL1, KIR2DS1, KIR2DS3, KIR2DS5",
        "HGNC:6329, HGNC:6333, HGNC:6335, HGNC:6337", "", "", CD158,

    "CD45RA", "HI100", "", "PTPRC", "PR_000001015", "", "Isoform of CD45",
    "CD45RO", "UCHL1", "", "PTPRC", "PR_000001017", "", "Isoform of CD45",

    "Ig light chain lambda", "MHL-38", "IGL", "HGNC:5853", "PR:000050185", "",
        "",

    "IgG Fc", "M1310G05", "IGHG1,IGHG3", "HGNC:5525, HGNC:5527", "", "", IgG,

    "c-Met", "12.1", "MET", "HGNC:7029", "", "", "",

    "CLA", "HECA-452", "SELPLG", "HGNC:10722", "", "", ""
)
