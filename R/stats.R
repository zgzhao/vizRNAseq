#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: stat_enrich.R
# Description:
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-03-07 16:38:44

#' The analysis requires two data sets: a gene list and a gene-to-ko mapper. A gene-to-ko list can be retrieved from KEGG website or obtained by sequence blast.
#'
#' KEGG pathway enrichment analysis.
#' @title KEGG pathway enrichment analysis.
#' @param gids character vector: gene ids to analysis
#' @param gmap data.frame or tibble with gene_id and ko (KEGG identifiers for genes such as K00185).
#' @param p.cut cut off value of BH adjust p value. Hypergeometric test is applied.
#' @param refGenes character vector of reference gene set.
#' @return tibble
#' @author ZG Zhao
#' @export
statKEGG <- function(gids, gmap, p.cut=0.05, refGenes=NULL) {
    stats <- KEGG_full_map(gmap)
    if(! is.null(refGenes))
        stats <- stats %>% filter(gene %in% {{refGenes}})
    
    stats <- stats %>% group_by(ko) %>%
        summarise(g.path=list(unique(gene))) %>%
        mutate(g.test=map(g.path, function(x) intersect(x, {{gids}}))) %>%
        mutate(n.path=map_int(g.path, length)) %>%
        mutate(n.test=map_int(g.test, length)) %>%
        filter(n.test > 0)

    allGenes <- unique(unlist(stats$g.path))
    lstGenes <- unique(unlist(stats$g.test))
    k <- length(lstGenes)   ## number of test
    nn <- length(allGenes)  ## total number of genes
    stats <- stats %>%
        mutate(p=pmap_dbl(list(q=n.test, m=n.path),
                          function(q, m) phyper(q, m, {{nn}} - q, {{k}}, lower.tail = FALSE)))
    stats$p.adj <- p.adjust(stats$p, method = "BH")
    stats %>% 
        left_join(
            KEGGFULL %>%
            select(ko, level, description) %>%
            distinct, by='ko') %>% 
        rename(n=n.test, N=n.path, genes=g.test, desc=description) %>%
        select(ko, level, n, N, p, p.adj, desc, genes) %>%
        mutate(ko=sub('^([0-9]+)$', 'KO:\\1', ko)) %>% 
        filter(p.adj < {{p.cut}}) %>% 
        arrange(p.adj, ko) %>%
        rename(ID=ko)
}

## map user-provided gene_id ~ go mapping to GO tree
oboMapGene <- function(gmap) {
    gmap <- gmap %>%
        group_by(go) %>%
        summarise(genes=list(gene_id))
    ans <- left_join(GOFULL, gmap, by = 'go') 
    for (i in 1:nrow(ans)) {
        ggs <- unlist(ans$genes[[i]])
        if(! is.null(ggs)) {
            gos <- na.omit(ans$isa[[i]])
            for(j in which(ans$go %in% gos)) {
                ggx <- unlist(ans$genes[[j]])
                ans$genes[[j]] <- union(ggx, ggs)
            }
        }
    }
    ans %>% mutate(gn=map_int(genes, length)) %>%
        filter(gn > 0) %>%
        select(-gn)
}

setGenes <- function(x, gset) {
    if(is.null(gset)) x else intersect(x, gset)
}

#' Gene ontology enrichment analysis
#'
#' Gene ontology enrichment analysis
#' @title GO terms enrichment/statistical analysis.
#' @param gids character vector of gene ids.
#' @param gmap data.frame or tibble with at least two columns: gene_id and go (GO identifiers).
#' @param refGenes character vector. Reference gene set, if NULL (default), all genes in the gene map (provided by gmap param) will be used.
#' @param p.cut numeric (0, 1). Default: 0.05
#' @return tibble, statistical results
#' @author ZG Zhao
#' @export
statGO <- function(gids, gmap=ATHGENEGO, refGenes=NULL, p.cut=0.05) {
    cat('Calculating, please wait a minute ... ...\n')
    if(! is.null(refGenes)) {
        gids <- intersect(gids, refGenes)
        gmap <- gmap %>%
            filter(gene_id %in% refGenes)
    }
    obo <- oboMapGene(gmap)
    galls <- unique(unlist(obo$genes))
    gtest <- intersect(galls, gids)
    
    res1 <- obo %>%
        mutate(genes = map(genes, setGenes, gset={{gids}})) %>% 
        mutate(n.test = map_int(genes, length))
    res2 <- obo %>% 
        mutate(genes = map(genes, setGenes, gset={{refGenes}})) %>% 
        mutate(n.ref = map_int(genes, length)) %>% 
        select(go, n.ref)
    nn <- length(galls)
    k <- length(gtest)
    stats <- left_join(res1, res2, by='go') %>%
        filter(n.ref > 0) %>%
        mutate(p=pmap_dbl(list(q=n.test, m=n.ref),
                          function(q, m) phyper(q, m, {{nn}} - q, {{k}},
                                                lower.tail = FALSE)))
    stats$p.adj <- p.adjust(stats$p, method = "BH")
    stats %>% filter(p.adj < p.cut) %>%
        select(space, go, n.test, n.ref, p, p.adj, desc, genes, everything(), -obsolete) %>% 
        rename(n=n.test, N=n.ref) %>% 
        arrange(p.adj) %>%
        rename(ID=go)
}
