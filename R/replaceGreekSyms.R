# replaceGreekSyms ----
#'@title replaceGreekSyms
#'@description Replace Greek symbols with either single letter or word.  Note
#'that this only replaces greek symbols I have encountered in antibody names to
#'avoid introducing errors by e.g. replacing "mu" in "immunoglobin".
#'@param x Character vector (n) containing Greek symbols to be replaced.
#'@param replacement Either "sym2letter" (replace Greek symbols with single
#'lowercase letter) or "sym2word" (replace Greek symbols with lowercase word,
#'e.g. "alpha") or "word2letter" (replace lowercase Greek symbol names with
#'single lowercase letters e.g. "alpha" to "a")
#'@return Character vector (n) with Greek symbols or names replaced
#'@importFrom stringr str_replace_all
#'@export
replaceGreekSyms <- function(x,
                        replace = c("sym2letter", "sym2word", "word2letter")){
    # Select replacement vector
    replace <- match.arg(replace)
    replace_to_vec = list("sym2letter" = GREEKSYM_TO_LETTER,
                       "sym2word" = GREEKSYM_TO_NAME,
                       "word2letter" = GREEKNAME_TO_LETTER)
    replace_vec <- replace_to_vec[[replace]]
    stringr::str_replace_all(x, replace_vec, names(replace_vec))
}
