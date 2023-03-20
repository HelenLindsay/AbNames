# Strip commonly observed prefixes and suffixes from ADT names
#
# A function that works with most of the data sets we've processed
#
# x = a vector of names
# anti = remove prefix "anti"?  Default: TRUE
# split_frac = fraction of entries that should include a potential delimiter in
# order to split entries.  Entries are often name_clone or name_gene,
# so isotype controls may have a different format
# Returns a data.frame of original and preprocessed names
# TO DO: proper loading of data set
#'@importFrom dplyr across
#'@importFrom dplyr matches
preprocessNames <- function(x, anti = TRUE, split_frac = 0.8){
    prefix1 = paste("ADT", "PROT", sep = "|")
    prefix2 = ifelse(isTRUE(anti), "[Aa]nti", "")
    suffix = paste("Ab", "PROT", "[ACTG]+", "TotalSeq[ABCD]",
                   "TS[ABCD]", sep = "|")
    prefix <- sprintf("^((%s)[[:punct:]]+)?(%s)?", prefix1, prefix2)
    suffix <- sprintf("([[:punct:]]+(%s))?$", suffix)

    y <- gsub(suffix, "", gsub(prefix, "", x))

    # Split if most entries of x have 2+ pieces using the same delimiter
    punct <- stringr::str_extract_all(y, "[[:punct:]]")
    frac_w_delim <- table(unlist(lapply(punct, unique))) / length(y)
    poss_delim <- names(frac_w_delim)[frac_w_delim > split_frac]
    result <- data.frame(Original = x, value = y)
    if (length(poss_delim) > 1){
        result <- result %>%
            dplyr::mutate(value_sp = strsplit(.data$value, poss_delim),
                          n_split = lengths(.data$value_sp),
                          # HERE CALC MODE SPLIT
                          value_sp = lapply(.data$value_sp, unique)) %>%
            tidyr::unnest_wider(.data$value_sp, names_sep = "_")

        # Quick match in the totalseq data set
        data(totalseq,  envir = environment())
        totalseq <- totalseq %>%
            dplyr::select(dplyr::all_of("Antigen", "Clone")) %>%
            unique()

        result <- result %>%
            # Note that totalseq doesn't have NA in Antigen or Clone cols
            dplyr::mutate(across(matches("value_sp"),
                ~ifelse(toupper(.x) %in% toupper(totalseq$Antigen),
                        .x, NA_character_))) %>%
            dplyr::mutate(value_sp = coalesce(across(matches("value_sp"))))
        # HERE IF_ELSE unusual split or value_sp is na -
        # value_sp not na - value_sp



         #                        .names = "{.col}_is_antigen"),
         #          across(matches("value_sp_[0-9]$"),
         #                 ~.x %in% totalseq$Clone,
         #                 .names = "{.col}_is_clone"))

        result_summary <- result %>%
            dplyr::summarise(across(matches("_is_"), sum)) / nrow(result)

        # If one of the columns frequently matched an antigen,
        # assume that this column is the antigen unless the split was unusual
        if (max(result_summary %>%
                dplyr::select(matches("antigen"))) > split_frac){
            # Cases: more than typical number of splits:

            #

        }


    }

    # Check for same antigen used multiple times - e.g. CD3.1
    result$temp_val <- .gsubNA("[[:punct:]][1-9]$", "", y)
    result$temp_num <- .gsubNA(".*[[:punct:]]([1-9]$)", "\\1", y)

    # TO DO: should this be temp_val?
    result$dups <- .dups(temp) & ! is.na(temp)
    # Because lots of gene names also end with numbers, we will require all
    # sequential numbers to be present
    # (This will not detect entries that have been split, but it becomes
    # complicated to figure out whether the entry is an antibody or clone)
    result <- result %>%
        dplyr::group_by(dplyr::all_of("temp_val")) %>%
        dplyr::mutate(n_exp = sum(.data$dups) == dplyr::n() & dplyr::n() > 1,
                      num_exp =
                          ifelse(is.na(.data$temp_num), NA,
                                 sort(.data$temp_num) == 1:dplyr::n()),
                      value = ifelse(.data$n_exp & .data$num_exp,
                                     .data$temp_val, .data$value)) %>%
        dplyr::select(dplyr::all_of("Original", "value"))

    return(result)

}
