fillByGroup_data <- function(){
    df <- data.frame(A = c(rep("A", 3), NA, rep("B", 4), rep("C", 3), NA),
                     B = c(rep("A", 2), NA, "A", rep("B", 4), rep("C", 3), NA),
                     C = c(NA, 1, 1, 1, 2, 3, NA, 2, 4, 4, NA, 5))

    df_c <- df %>% dplyr::filter(! is.na(A) & ! is.na(B))
    df_d <- c(6, NA, 6, 7, 8, NA, 8, NA, 9, 8, NA, 10)

    # Simple filling function using tidyr::fill
    partial_f <- function(df, col, gp){
        df %>%
            dplyr::group_by(!!!syms(gp)) %>%
            tidyr::fill(!!sym(col), .direction = "updown") %>%
            dplyr::ungroup()
    }

    return(list(df = df, df_c = df_c, df_d = df_d, partial_f = partial_f))
}
