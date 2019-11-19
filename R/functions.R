#' Nuclear genes
#'
#' Get arabidopsis nuclear genes from a vecter.
#' @title nugenes function
#' @param x character vecter.
#' @return character vecter.
#' @author ZG Zhao
#' @export
nugenes <- function(x) {
    grep("^AT[1-5]G", x, value = TRUE)
}

#' Multiple intersect
#'
#' Get intersect of a list. The atoms exist in all list elements are returned.
#' @title Multiple intersect function
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
