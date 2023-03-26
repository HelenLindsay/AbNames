# Furnctions returning vectors for replacing Greek symbols or names ----

GREEKSYM_TO_NAME <- function(){
        return(c("\u03b1" = "alpha", "\u03b2" = "beta", "\u03b3" = "gamma",
                 "\u03b4" = "delta", "\u03b5" = "epsilon",
                 "\u03b6" = "zeta", "\u03ba" = "kappa", "\u03bb" = "lambda",
                 "\u03c4" = "tau", "\u03bc" = "mu"))
}

GREEKSYM_TO_LETTER <- function(){
        return(c("\u03b1" = "a", "\u03b2" = "b", "\u03b3" = "g",
                 "\u03b4" = "d", "\u03b5" = "e", "\u03b6" = "z",
                 "\u03ba" = "k", "\u03bb" = "l", "\u03c4" = "t",
                 "\u03bc" = "m"))
}

# Mu is not included here because it is in e.g. immunoglobin
GREEKNAME_TO_LETTER <- function(){
        c(alpha = "a", beta = "b", gamma = "g", delta = "d",
          epsilon = "e", zeta = "z", kappa = "k", lambda = "l",
          Alpha = "a", Beta = "b", Gamma = "g", Delta = "d",
          Epsilon = "e", zeta = "z", Kappa = "k", Lambda = "l")
}


# Return a table of known duplicated antibody clone names ----
CLONE_DUPS <- function(){
    vendors <- c("BioLegend", "BD Biosciences")

    df <- data.frame(Clone = c("3D12", "3D12", "3D12", "3D12", "3D12"),
                     Antibody = c("HLA-E", "HLA-E", "HLA-E", "HLA-E",
                                  "CCR7 (CD197)"),
                     Vendor = vendors[1, 1, 1, 1, 2],
                     TotalSeq_Cat = c("A", "B", "C", "D", NA),
                     Cat_Number = c("342617", "342621", "342619", "342623",
                                    "940014"))
}

# Return a table of isotype controls (checked manually)
ISOTYPE_CONTROLS <- function(){
    return(data.frame(Antigen = c("IgG2a", "IgG2a", "IgG1",
                                  "IgG", "IgG1", "IgG2b",
                                  "IgG1", "IgG2c", "IgG2b"),
                      Clone = c("MOPC-173", "RTK2758", "G0114F7",
                                "HTK888", "MOPC-21", "MPC-11",
                                "RTK2071", "RTK4174", "RTK4530")))
}
