#' Construct KEGG metabolic network from a given starting point.
#'
#' Giving an KEGG identifier
#' @title Induce KEGG network with a given root node.
#' @param root character string, KEGG identifier from which the KEGG network will be subset.
#' @return tbl_graph 
#' @author ZG Zhao
#' @export
induce_KEGG_network <- function(root) {
    dx <- filter(KEGGFULL, ko %in% {{root}}) 
    if (nrow(dx) < 0)
        stop('Not a valid KEGG identifier!')
    kk <- dx$ko[1]
    ll <- dx$level[1]
    csels <- LETTERS[1:4]
    csels <- csels[csels >= ll]
     
    knodes <- KEGGFULL %>% 
        filter(ko == {{kk}} | .data[[ll]]=={{kk}}) %>% 
        rename(D=ko) %>% select(all_of(csels))
    ans <- NULL
    for (i in 2:ncol(knodes)) {
        anx <- knodes[, (i-1):i]
        colnames(anx) <- c('from', 'to')
        ans <- rbind(ans, anx)
    }
    ans %>% drop_na %>% 
        distinct(from, to) %>% 
        as_tbl_graph %>% 
        activate('nodes') %>% 
        left_join(KEGGFULL %>% 
                      rename(name=ko) %>% 
                      distinct(name, .keep_all=TRUE),
                  by='name') %>% 
        select(name, level, description)
}
