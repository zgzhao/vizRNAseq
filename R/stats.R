#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: stat_enrich.R
# Description:
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-03-07 16:38:44

.statOnto <- function(gids, gmap, allGids, refGids=allGids, ngene.min=10) {
    ## NOTE: the input gmap has two columns: ID and genes
    ans <- gmap %>% 
        mutate(
            g.path=map(genes, function(x) intersect(x, {{allGids}})),
            g.test=map(genes, function(x) intersect(x, {{gids}})),
            g.refs=map(genes, function(x) intersect(x, {{refGids}})),
            n.path=map_int(g.path, length),
            n.test=map_int(g.test, length),
            n.refs=map_int(g.refs, length)) %>%
        filter(n.path >= {{ngene.min}}, n.test > 0)
    
    pathGenes <- unique(unlist(ans$g.path))
    testGenes <- unique(unlist(ans$g.test))
    k <- length(testGenes)  ## number of test genes
    n <- length(pathGenes)  ## number of all genes
    ## hypergeometric test
    ans$p <- apply(ans[, c('n.test', 'n.path')],
                1, FUN=function(x) {
                    q <- x[1]
                    m <- x[2]
                    phyper(q-1, m, n - m, k, lower.tail = FALSE)
                })
    ans %>%
        select(ID, n.test, n.path, n.refs, p, g.test) %>% 
        rename(n=n.test, N=n.path, M=n.refs, genes=g.test)
}

#' The analysis requires two data sets: a gene list and a gene-to-ko mapper. A gene-to-ko list can be retrieved from KEGG website or obtained by sequence blast.
#'
#' KEGG pathway enrichment analysis.
#' @title KEGG pathway enrichment analysis
#' @param gids character vector: gene ids to analysis
#' @param gmap data.frame or tibble with gene_id and ko (KEGG identifiers for genes such as K00185).
#' @param p.cut cut off value of BH adjust p value. Hypergeometric test is applied.
#' @param allGids character vector of reference gene set.
#' @param ngene.min integer, minimum number of genes in test path (default: 10).
#' @return tibble
#' @author ZG Zhao
#' @export
statKEGG <- function(gids, gmap, allGids, refGids=allGids, p.cut=0.05,
                     ngene.min=10, plant=TRUE) {
    ## NOTE: do all filters before statistical test!
    gids <- intersect(gids, allGids)
    gmap <- onto_full_map(gmap, onto.type='ko')
    if (plant) {
        ko_common <- c('09100', '09120', '09130', '09140', '09150')
        kos <- KEGGFULL %>%
            mutate(sels=map_lgl(ancestors, function(x) length(intersect(x, {{ko_common}})) > 0)) %>% 
            filter(ID %in% {{ko_common}} | sels) %>%
            pull(ID)
        gmap <- filter(gmap, ID %in% kos)
    }
    ## hypergeometric test api
    ans <- .statOnto(gids, gmap, allGids, refGids, ngene.min) 
    ans$p.adj <- p.adjust(ans$p, method = "BH")
    ans %>% 
        filter(p.adj < {{p.cut}}) %>% 
        left_join(
            KEGGFULL %>%
            select(ID, level, desc) %>%
            distinct, by='ID') %>% 
        select(ID, level, n, N, M, p, p.adj, desc, genes) %>%
        arrange(p.adj, ID)
}

#' Gene ontology enrichment analysis
#'
#' Gene ontology enrichment analysis
#' @title GO terms enrichment/statistical analysis
#' @param gids character vector of gene ids.
#' @param gmap data.frame or tibble with at least two columns: gene_id and go (GO identifiers).
#' @param allGids character vector. Reference gene set, if NULL (default), all genes in the gene map (provided by gmap param) will be used.
#' @param p.cut numeric (0, 1). Default: 0.05
#' @return tibble, statistical results
#' @author ZG Zhao
#' @export
statGO <- function(gids, gmap, allGids, refGids, p.cut=0.05, ngene.min=10) {
    cat('Calculating, please wait a minute ... ...\n')

    ## NOTE: do all filters before statistical test!
    gids <- intersect(gids, allGids)
    gmap <- onto_full_map(gmap, onto.type='go')
    if(missing(refGids)) refGids <- allGids

    ## hypergeometric test api
    ans <- .statOnto(gids, gmap, allGids, refGids, ngene.min)
    ans$p.adj <- p.adjust(ans$p, method = "BH")

    ans %>% filter(p.adj < {{p.cut}}) %>%
        left_join(GOFULL, by='ID') %>% 
        arrange(space, p.adj, ID) %>%
        select(ID, space, n, N, M, p, p.adj, desc, genes)
}


#' Convert enrichment results to graph/network data.
#'
#' Refer to \code{\link[tidygraph:as_tbl_graph]{tidygraph::as_tbl_graph}} for more info about `tbl_graph`
#' @title network construction with enrichment data
#' @param df_stat tibble, results from \code{\link{statKEGG}} or \code{\link{statGO}}
#' @param group_label character, set this if you want to join graph later with \code{\link[tidygraph:graph_join]{tidygraph::graph_join}}
#' @return tbl_graph object
#' @author ZG Zhao
#' @export
enrich_to_graph <- function(df_stat, group_label='a') {
    if(nrow(df_stat) < 1)
        return(tibble())
    fullmap <- if(grepl('^GO:', df_stat$ID[1])) GOFULL else KEGGFULL
    ## create edges
    xedges <- df_stat %>% select(ID) %>% 
        left_join(fullmap %>% select(ID, isa), by = 'ID') %>%
        select(isa, ID) %>%
        rename(from=isa, to=ID) %>%
        unnest(from) %>%
        filter(from %in% df_stat$ID) %>%
        left_join(
            df_stat %>% select(ID, genes) %>%
            rename(from=ID, g1=genes),
            by = 'from') %>% 
        left_join(
            df_stat %>% select(ID, genes) %>%
            rename(to=ID, g2=genes),
            by = 'to') %>%
        mutate(x.shared=pmap_dbl(list(g1, g2), function(a, b) {
            length(unique(b))/length(unique(a))
        }))
    ## set ratio shared genes (by parents)
    xlen <- function(x) length(unique(unlist(x)))
    xedges <- xedges %>%
        group_by(from) %>% 
        summarise(p.shared = xlen(g2)/xlen(g1)) %>%
        right_join(xedges, by = 'from', multiple = 'all') %>%
        select(from, to, x.shared, p.shared)
    ## create graph now
    xgraph <- xedges %>% as_tbl_graph
    xnodes <- setdiff(df_stat$ID, unlist(xedges[, 1:2]))
    if (length(xnodes) > 0)
        xgraph <- xgraph %>% bind_nodes(tibble(name=xnodes))
    ## set node data
    pshared <- xedges %>%
        select(to, p.shared) %>%
        rename(ID=to) %>%
        group_by(ID) %>%
        summarise(p.shared=max(c(0, p.shared), na.rm=TRUE))
    df_stat <- df_stat %>%
        mutate(m.shared=pmap_dbl(list(n, M), function(a, b) a/b)) %>% 
        left_join(pshared, by = 'ID') %>%
        mutate(p.shared=ifelse(is.na(p.shared), 0, p.shared)) %>% 
        rename(name=ID)
    ## bind node data
    xgraph %>% activate('nodes') %>%
        left_join(df_stat %>% select(-genes), by='name') %>%
        mutate(group={{group_label}},
               label=name,
               name=paste(label, group, sep='_'))
}


