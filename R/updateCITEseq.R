# Add new data to the CITEseq data set and fill missing information if possible
# To do: should we keep any extra columns not in citeseq?
# x: data to add to CITEseq data set
# extra: logical(1) should columns not in CITEseq be kept?
# name: study name
#'@importFrom rlang check_installed
updateCITEseq <- function(x, name = NULL, extra = FALSE){
    if (! rlang::check_installed(stringi)){
        stop("This function requires stringi to be installed")
    }

    citeseq_fname <- system.file("extdata", "citeseq.csv", package = "AbNames")
    citeseq <- read.csv(citeseq_fname)

    # Study name must be provided if not already in x ----
    if (! "Study" %in% colnames(x)){
        if (is.null(name)) {
            stop("Please provide a name for the new data")
        }
        x$Study <- name
    }

    # Antigen and Vendor columns must be present ----
    if (any(! c("Antigen", "Vendor") %in% colnames(x))){
        stop("Data to add must contain columns 'Antigen' and 'Vendor'")
    }

    # Should columns not in citeseq be kept? ----
    if (isFALSE(extra)){
        # Only keep columns already present in CITE-seq collection
        x <- unique(x[, intersect(colnames(x), colnames(citeseq))])
    }

    # (note cite-seq has greek syms...) -----------------
    # Check for non-ascii characters
    x$Antigen <- replaceGreekSyms(x$Antigen)
    if (any(grepl("[^ -~]", x$Antigen))){
        warning("Non-ascii characters replaced in Antigen column")
        x$Antigen <- stringi::stri_trans_general(x$Antigen, "latin-ascii")
    }
    # -------------------------------------------------------

    # Are there already ID columns present?
    # TO DO: NOT CURRENTLY IN CITEseq data!
    id_cols <- c("HGNC_ID", "ENSEMBL_ID", "HGNC_SYMBOL", "ENTREZ_ID")
    id_cols <- intersect(id_cols, colnames(x))

    # If ID columns are provided, check that entries are valid
    if (length(id_cols) > 0){

        # ---- to update --------
        data(hgnc)
        # -----------------------

        # Do this as a loop to know which is inconsistent?
        inconsistent <- dplyr::anti_join(x, hgnc, by = id_cols)
        if (nrow(inconsistent) > 0){
            stop("ID rows are inconsistent")
        }
        # Update with correct information
    }

    # Update suggested antigen

    # If vendor information is provided, fill missing and check consistency ----
    vendor_cols <- c("Cat_Number", "Vendor", "Oligo_ID", "Clone",
                     "RRID", "TotalSeq_Cat")
    vendor_cols <- Reduce(intersect, list(vendor_cols,
                                          colnames(x),
                                          colnames(citeseq)))

    # Format x vendor cols to match citeseq
    x <- x %>%
        dplyr::mutate(across(all_of(vendor_cols), as.character)) %>%
        # remove unnecessary whitespace
        dplyr::mutate(across(where(is.character), stringr::str_squish))

    # Vendor must match, does x have multiple vendors?

    # Select complete cases from citeseq data where any column matches x
    #cs <- citeseq %>%
    #    tidyr::drop_na(dplyr::any_of(vendor_cols)) %>%
    #    union_join(x %>% dplyr::select(vendor_cols))

}




checkCITEseq <- function(citeseq){
    # One Cat_Number should match one Oligo / Vendor / ID / Clone / RRID
    # (but check e.g. Cat_Number is "custom made")

    # 1 Cat_Number -> 1 combination of Vendor, Oligo_ID, Clone
    # 1 RRID -> 1 combination of Vendor, Oligo_ID, Clone
    #

    # Same suggested Antigen / Clone / Vendor -> One Gene Id?

    #dplyr::group_by(Suggested_Antigen, Oligo_ID, Clone, TotalSeq_Cat) %>%
    #tidyr::fill(c(Cat_Number, Reactivity), .direction = "updown") %>%

    x <- citeseq %>%
        dplyr::group_by(Cat_Number) %>%
        dplyr::mutate(nrrid = n_distinct(RRID, na.rm = TRUE),
                      nclone = n_distinct(tolower(Clone), na.rm = TRUE),
                      ntotalseq = n_distinct(TotalSeq_Cat, na.rm = TRUE),
                      noligo = n_distinct(Oligo_ID, na.rm = TRUE))



    y <- x %>% dplyr::filter(nrrid > 1 | nclone > 1 | ntotalseq > 1 |
                                 noligo > 1 | nbarcode > 1) %>%
        dplyr::pull(Cat_Number) %>% unique()

    z <- x %>% dplyr::filter(Cat_Number %in% y) %>%
        dplyr::select(Antigen, Study, Cat_Number, Clone, nclone, TotalSeq_Cat,
                      ntotalseq, Oligo_ID, noligo, Barcode_Sequence, nbarcode,
                      RRID, nrrid) %>%
        dplyr::filter(! is.na(Cat_Number) & ! Cat_Number == "custom made") %>%
        dplyr::arrange(Cat_Number)

}


fillCITEseq <- function(citeseq){
    # If RRID matches, fill Clone, Oligo_ID, Cat_Number, Vendor?


    # If Cat_Number matches, fill Clone, Oligo_ID, Cat_Number, Vendor?

    # If Vendor and Cat_Number match, fill Oligo_ID, Clone
    # Fill information from Cat_Number if available
    dplyr::group_by(Cat_Number) %>%
        tidyr::fill(c(Clone, Oligo_ID, TotalSeq_Cat, Reactivity),
                    .direction = "updown")


    # Select the most common Selected_Antigen, replace all members of group
    dplyr::add_count(Cat_Number, Suggested_Antigen) %>%
        dplyr::group_by(Cat_Number) %>%
        dplyr::mutate(n = ifelse(n == 1 & grepl("kappa", Suggested_Antigen),
                                 100, n)) %>% # Prefer kappa if it's there
        dplyr::mutate(SA = Suggested_Antigen[n == max(n)][1],
                      Suggested_Antigen = ifelse(is.na(Cat_Number),
                                                 Suggested_Antigen, SA),
                      Control = TRUE) %>%

        # Fill Cat_Number by grouping with Suggested_Antigen
        # (instead of Antigen as in later filling code)
        # (Note: have manually checked for NA groups)
        dplyr::group_by(Suggested_Antigen, Oligo_ID, Clone, TotalSeq_Cat) %>%
        tidyr::fill(c(Cat_Number, Reactivity), .direction = "updown") %>%
        dplyr::ungroup() %>%
        dplyr::select(-n, -SA)

    na_cat <- all_clones %>%
        dplyr::filter(is.na(Cat_Number))

    # Note: have checked that Clone is never NA if Cat_Number is "custom made"
    all_clones <- all_clones %>%
        dplyr::filter(! is.na(Cat_Number)) %>%
        dplyr::group_by(Cat_Number) %>%
        tidyr::fill(c(Clone, Oligo_ID, RRID), .direction = "updown") %>%
        dplyr::ungroup() %>%
        dplyr::full_join(na_cat)


    # Fill in missing catalogue numbers if Oligo_ID, Clone and Antigen match:
    # TotalSeq_Cat is needed, it is possible to share oligo, antigen and clone
    # but not cat number
    # Hao_2021 and Liu_2021 have entries for some but not all Cat_Numbers
    # It appears info can be added for Liu, but Hao entries are custom
    # (But do not ever have other group members)

    na_oligo_clone_ag <- all_clones %>%
        dplyr::filter(is.na(Oligo_ID) | is.na(Clone) |
                          is.na(Antigen) | is.na(TotalSeq_Cat) |
                          grepl("custom", Cat_Number))

    all_clones <- all_clones %>%
        dplyr::filter(! is.na(Oligo_ID) & ! is.na(Clone) &
                          ! is.na(Antigen) & ! is.na(TotalSeq_Cat) &
                          ! grepl("custom", Cat_Number)) %>%
        dplyr::group_by(Antigen, Oligo_ID, Clone, TotalSeq_Cat) %>%
        tidyr::fill(Cat_Number, .direction = "updown") %>%
        dplyr::ungroup() %>%
        dplyr::full_join(na_oligo_clone_ag)




}



# Check RRID and Cat_Number correct if info available

#citeseq %>%
# group_by(RRID) %>%
# dplyr::filter(n_distinct(HGNC_ID) > 1 |
# n_distinct(Clone) > 1 | n_distinct(Cat_Number) > 1)

# If symbols are provided, check that symbols are valid
# Fill citeseq RRID by Cat_Number


# Bind to CITE-seq

# If all rows have an identifier, group by ID and select majority

# If has catalogue number,
