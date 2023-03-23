.prettyLabel <- function(x, min.len=10, min.word=2) {
    x <- trimws(x)
    words <- strsplit(x, ' +')[[1]]
    nw <- length(words)
    ncs <- nchar(words)
    xmin <- nchar(x)
    if (xmin < min.len || nw < min.word)
        return(x)
    nn <- 0
    for (i in 1:(nw-1)) {
        ndff <- abs(sum(ncs[1:i]) - sum(ncs[(i+1):nw]))
        if(ndff < xmin) {
            nn <- i
            xmin <- ndff
        }
    }
    aa <- paste(words[1:nn], collapse=' ')
    bb <- paste(words[(nn+1):nw], collapse=' ')
    paste(aa, bb, sep='\n')
}

.level2size <- function(x) {
    ss <- c(30, 24, 18, 12)
    names(ss) <- LETTERS[1:4]
    ss[x]
}


getEdges <- function(g) {
    fun1 <- function(x) sub('\\|.+$', '', x)
    fun2 <- function(x) sub('^.+\\|', '', x)
    ans <- tibble(name = attr(E(g), 'vnames')) %>%
        mutate(from=map_chr(name, fun1),
               to=map_chr(name, fun2))
    atts <- edge_attr(g)
    if (length(atts) > 0)
        ans <- cbind(ans,  atts %>% as_tibble)
    ans
}

#' Plot network of created from KEGG/GO enrichment results
#'
#' Joined tbl_graphs are support. Please pipe `facet_nodes` after `plotEnrich`. 
#' @title Visualization of enrich network
#' @param egraph tbl_graph, results from function `enrich_graph`.
#' @param layout param passed to `layout` for `ggraph`
#' @param group.cut numeric (0~1), cutoff value to determine whether genes in children nodes are dominated in a parent.
#' @param path.cut  numeric (0~1), cutoff value to determine whether genes in a node are dominated in the same pathway.
#' @param x.expand numeric (0~1), for `expand` of `scale_x_continuous`.
#' @param y.expand numeric (0~1), for `expand` of `scale_x_continuous`.
#' @return ggraph
#' @author ZG Zhao
#' @export
plotEnrich <- function(egraph, layout='kk', 
                       group.cut=0.5, path.cut=0.8, label.desc=TRUE,
                       node.size.factor=10, node.size.cex=1, node.label.cex=1,
                       node.colors=c('gray50', 'red'),
                       edge.colors=c('gray70', 'orange'),
                       edge.width=1,
                       edge.arrow.angle=30, edge.arrow.length=3, edge.arrow.type='open', 
                       x.expand=0, y.expand=0) {
    
    gg <- egraph %>%
        activate('nodes') %>% 
        mutate(
            m.shared = n/M,
            size=log(n + {{node.size.factor}}, {{node.size.factor}}) * {{node.size.cex}},
            node_color=ifelse(m.shared>{{path.cut}}, 'AEGS', 'non-AEGS'), 
            node_color=factor(node_color, levels=c('non-AEGS', 'AEGS')),
            label=sub('_.+$', '', name)) %>% 
        activate('edges') %>%
        mutate(edge_color=ifelse(p.shared > group.cut, 'to CEGS', 'to non-CEGS'),
               edge_color=factor(edge_color, levels=c('to non-CEGS', 'to CEGS')))
    ## set labels
    if(label.desc) {
        gg <- gg %>%
            activate('nodes') %>% 
            mutate(
            label=pmap_chr(list(label, desc), function(kk, dd) {
                kk <- sub('_.+$', '', kk)
                dd <- sub("^.+;", '', dd)
                paste(kk, dd, sep=": ")
            }),
            label=map_chr(label, .prettyLabel))
    }

    gg %>% 
        ggraph(layout = layout) +
        geom_node_point(aes(color=node_color, size=I(size))) +
        geom_edge_link(aes(color=edge_color),
                       arrow=arrow(angle=edge.arrow.angle,
                                   length=unit(edge.arrow.length, 'mm'),
                                   type=edge.arrow.type),
                       width=edge.width) +
        geom_node_text(aes(label=label), size=5*{{node.label.cex}}) +
        scale_x_continuous(expand = c(x.expand, x.expand)) +
        scale_y_continuous(expand = c(y.expand, y.expand)) +
        scale_color_manual(values = alpha(node.colors, 0.5)) +
        scale_edge_color_manual(values = alpha(edge.colors, 0.7))
}



#' @export
basalGOnet <- function(enodes, category=c('generic', 'goslim_plant')) {
    xgraph <- GOFULL %>%
        mutate(cat=map_lgl(cat, function(x) {{category}} %in% x)) %>%
        filter(cat, !obsolete) %>%
        rename(from=isa, to=ID) %>%
        select(from, to) %>%
        unnest(from) %>%
        as_tbl_graph %>%
        mutate(core=name %in% {{enodes}})
    ne1 <- 10
    ne2 <- 0
    while(ne1 > ne2) {
        ne1 <- xgraph %>% 
            as_tibble %>% nrow
        xgraph <- xgraph %>%
            mutate(d=centrality_degree(mode = "out")) %>%
            filter(!(! core & d == 0))
        ne2 <- xgraph %>%
            as_tibble %>% nrow
    }
    xgraph %>%
        select(-d)
}
