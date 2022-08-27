# Strip commonly observed prefixes and suffixes from ADT names
#
# A function that works with most of the data sets we've processed
#
# x = a vector of names
# anti = remove prefix "anti"?  Default: TRUE
# split_frac = fraction of entries that should include a potential delimiter in
# order to split entries.  Entries are often name_clone or name_gene,
# so isotype controls may have a different format
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
    if (length(poss_delim) > 1){
        result <- data.frame(Original = x, value = strsplit(y, poss_delim)) %>%
            tidyr::unnest(value)
    } else {
        result <- data.frame(Original = x, value = y)
    }

    # Check for same antigen used multiple times - e.g. CD3.1
    result$temp_val <- .gsubNA("[[:punct:]][1-9]$", "", y)
    result$temp_num <- .gsubNA(".*[[:punct:]]([1-9]$)", "\\1", y)
    result$dups <- .dups(temp) & ! is.na(temp)
    # Because lots of gene names also end with numbers, we will require all
    # sequential numbers to be present
    # (This will not detect entries that have been split, but it becomes
    # complicated to figure out whether the entry is an antibody or clone)
    result <- result %>%
        dplyr::group_by(temp_val) %>%
        dplyr::mutate(n_exp = sum(dups) == dplyr::n() & dplyr::n() > 1,
                      num_exp = ifelse(is.na(temp_num), NA,
                                             sort(temp_num) == 1:dplyr::n()),
                      value = ifelse(n_exp & num_exp, temp_val, value)) %>%
        dplyr::select(Original, value)

}
