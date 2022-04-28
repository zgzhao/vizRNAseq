
#' Calculate shared genes of each pair of GO terms or KEGG pathways.
#'
#' Calculate shared genes of each pair of GO terms or KEGG pathways.
#' @title Calculate shared genes.
#' @param df data.frame/tibble returned by `statGO` or `statKEGG`
#' @param percent TRUE/FALSE. If TRUE, return percentage of genes (union of two nodes).
#' @return tibble with three columns: v1, v2 and w (number or percent of shared genes).
#' @author ZG Zhao
#' @export
sharedGenes <- function(df, percent=FALSE) {
    genes <- df$genes
    nn <- nrow(df)
    ans <- matrix(0, ncol = nn, nrow = nn)
    for(i in 1:nrow(df)) {
        gg1 <- genes[[i]]
        for(j in 1:nrow(df)) {
            if(j > i) next
            gg2 <- genes[[j]]
            nxx <- length(intersect(gg1, gg2))
            if(percent) nxx <- nxx/length(union(gg1, gg2))
            ans[i, j] <- nxx
        }
    }
    vss <- df$ID
    colnames(ans) <- vss
    ans <- as_tibble(ans) %>%
        mutate(v1 = {{vss}}) %>%
        pivot_longer(-v1, names_to = 'v2', values_to = 'w') %>%
        filter(w > 0, v1 != v2)
    ans
}

#' Calculate shared genes and create a comunity network.
#'
#' Use \code{\link{sharedGenes}} for a table of shared genes only.
#' @title Create shared gene network/graph.
#' @param df data.frame/tibble returned by `statGO` or `statKEGG`
#' @param percent TRUE/FALSE. Use percentage or number of shared genes as edge weights.
#' @return list. Element 'graph': igraph object. Element 'community': igraph group object.
#' @author ZG Zhao
#' @export
graphTerms <- function(df, percent=TRUE) {
    ans <- sharedGenes(df, percent)
    gg <- graph.data.frame(ans, directed = FALSE)
    kc <- fastgreedy.community(gg, weights = E(gg)$w)
    list(graph = gg, community= kc)
}

