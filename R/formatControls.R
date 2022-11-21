formatControls <- function(){

    control_manual <- tibble::tribble(~Antigen, ~Suggested_Antigen,
                                      "Arm.ham.IgG", "Armenian Hamster IgG",
                                      "IgG Fc", "IgG (Fc)",
                                      "IgG.Fc", "IgG (Fc)",
                                      "IgG2A", "IgG2a")


    rename_sp <- structure(c("Mouse ", "Rat ", "Rat ", "Mouse ", "Rat ", "kappa",
                             "lambda", "Armenian Hamster", "Hamster IgG"),
                           names = c("^m", "^r", "Rag_", "Mouse[_\\.]", "Rat[_\\.]",
                                     "[κk]$", "λ",  "Arm\\.ham\\.", "HamsterIgG"))

    controls <- all_clones %>%
        dplyr::filter(grepl("Mouse|Ra[tg]|[Hh]am|[Cc]ontrol|Ctrl|[Ii]so",
                            Antigen)) %>%

        # Standardise the Antigen name in Suggested_Antigen column
        dplyr::mutate(Suggested_Antigen =
                          gsub("^control[ _]([iI]sotype )?", "", Antigen),
                      Suggested_Antigen =
                          gsub(" [iI]sotype Ctrl$", "", Suggested_Antigen),
                      Suggested_Antigen =
                          gsub("^Iso|^Isotype control,? ", "", Suggested_Antigen),
                      Suggested_Antigen = gsub("\\.[Ii]so$|", "", Suggested_Antigen),
                      Suggested_Antigen = gsub("_[12]$", "", Suggested_Antigen),
                      Suggested_Antigen = stringr::str_replace_all(
                          Suggested_Antigen, rename_sp, names(rename_sp)))

}
