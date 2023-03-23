#' Filter nuclear genes from a gene name vector (remove plastid genes).
#'
#' This function is for Arabidopsis only.
#' @title Arabidopsis: filter nuclear genes
#' @param x character vecter.
#' @return character vecter.
#' @author ZG Zhao
#' @export
ath_nugenes <- function(x) {
    grep("^AT[1-5]G", x, value = TRUE)
}

#' Get intersect of a list. The atoms exist in all list elements are returned.
#'
#' The general \code{\link{intersect }} function accepts two vectors only. `mintersect` is more convenient for multiple intersect.
#' @title intersect function for list
#' @param x list.
#' @return vector
#' @author ZG Zhao
#' @export
mintersect <- function(x) {
    if (! is.list(x) ) stop("Input data must be a list.")
    res <- unlist(x)
    for (aa in x) {
        res <- intersect(res, aa)
    }
    res
}
