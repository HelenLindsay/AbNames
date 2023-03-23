test_that("annotateAb works", {
    expect_error(annotateAb(x = data.frame(Antibody = "CD158b")))

    df <- data.frame(Antigen = c("CD279",
                                 "CD279 (PD-1)",
                                 "CD279_PD1",
                                 "CD279(PD-1)",
                                 "PD-1 (CD279)",
                                 "PD1 (CD279)",
                                 "KLRG1 (MAFA)",
                                 "CD11a.CD18"))
    res <- annotateAb(df)
    same_ab <- c("CD279", "CD279 (PD-1)", "CD279_PD1", "CD279(PD-1)",
                  "PD-1 (CD279)", "PD1 (CD279)")
    same_ab_res <- res[match(same_ab, res$Antigen), "ALT_ID"]

    # Matching genes map to the same ID
    expect_true(length(unique(same_ab_res)) == 1)

    # Expect that there are multiple gene IDs for ab targeting multiple genes
    multi_gene <- c("KLRG1 (MAFA)", "CD11a.CD18")
    multi_gene_id <- res[match(multi_gene, res$Antigen), "ALT_ID"]
    expect_true(all(grepl(",", multi_gene_id)))
})
