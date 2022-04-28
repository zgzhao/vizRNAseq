#' Gene ontology network. 
#'
#' Create graph/network from gene ontology file (obo).
#' @title Create gene ontology graph/network 
#' @param obofile Character. obo file or one of 'full', 'generic', 'plant', 'pir', metagenomics'.
#' @return graphNEL object
#' @author ZG Zhao
#' @export
GO_graph <- function(obo = 'full') {
    if(file.exists(obo)) oinfo <- parseOBO(obo)
    else oinfo <- GOFULL
    obo2graph(oinfo)
}

#' @export
obo2graph <- function(oinfo) {
    ess <- oinfo %>%
        select(go, isa) %>% 
        unnest(isa) %>%
        drop_na(isa) %>%
        distinct %>%
        group_by(isa) %>%
        summarise(ess=list(unique(go)))
    exx <- pull(ess)
    names(exx) <- ess$isa
    nss <- unique(unlist(ess))
    gg <- graphNEL(nodes=nss, edgeL = exx, edgemode = 'directed')
    oinfo <- filter(oinfo, go %in% nodes(gg))
    nodeDataDefaults(gg, 'label') <- ''
    nodeData(gg, oinfo$go, 'label') <- oinfo$desc
    gg
}
