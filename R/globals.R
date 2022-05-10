# Vectors for replacing Greek symbols or names ----

# Rho?

GREEKSYM_TO_NAME = c("\u03b1" = "alpha", "\u03b2" = "beta", "\u03b4" = "gamma",
                     "\u03b4" = "delta", "\u03b5" = "epsilon",
                     "\u03b6" = "zeta", "\u03ba" = "kappa",
                     "\u03bb" = "lambda", "\u03c4" = "tau", "\u03bc" = "mu")

GREEKSYM_TO_LETTER =  c("\u03b1" = "a", "\u03b2" = "b", "\u03b4" = "g",
                        "\u03b4" = "d", "\u03b5" = "e", "\u03b6" = "z",
                        "\u03ba" = "k", "\u03bb" = "l", "\u03c4" = "t",
                        "\u03bc" = "m")

# Tau is not included, is it necessary?
# Mu is not included here because it is in e.g. immunoglobin
GREEKNAME_TO_LETTER = c(alpha = "a", beta = "b", gamma = "g", delta = "d",
                        epsilon = "e", zeta = "z", kappa = "k", lambda = "l",
                        Alpha = "a", Beta = "b", Gamma = "g", Delta = "d",
                        Epsilon = "e", zeta = "z", Kappa = "k", Lambda = "l")
