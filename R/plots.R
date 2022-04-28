
#' Construct KEGG network (graphNEL object) from a list of KOs.
#'
#' Duplicated and cyclic edges will be removed in the final graph.
#' @title Function kinfo2graph
#' @param kos character vector of KEGG indentifiers (KOs).
#' @param include.genes TRUE/FALSE. If true, k genes (level D entries in KEGG pathway) will be included.
#' @return an igraph object
#' @author ZG Zhao
#' @export
kinfo2graph <- function(kos, include.genes=FALSE){
    if(! include.genes) kos <- grep("^[0-9]", kos, value=TRUE)
    kos <- intersect(kos, keggs$ko)
    kegx <- keggs[keggs$ko %in% kos, c("A", "B", "C", "ko")]
    eds <- apply(kegx, 1, FUN=function(x){
        vv <- setdiff(x, NA)
        nn <- length(vv)
        if(nn < 2) return(NA)
        paste(vv[1:(nn - 1)], vv[2:nn])
    })
    eds <- unique(unlist(eds))
    eds <- setdiff(eds, NA)
    eds <- sapply(eds, FUN=function(x)strsplit(x, " ")[[1]])
    eds <- t(eds)
    vss <- unique(unlist(eds))
    g <- graphNEL(vss, edgemode = "directed")
    for(i in 1:nrow(eds)) g <- addEdge(eds[i, 1], eds[i, 2], g)
    g
}


#' Plot KEGG networks, mostly for results obtained by \code{\link{KEGG_enrich}}.
#'
#' Refer to \code{\link{graph::plot.graphNEL}} for more details.
#' @title Plot KEGG network or igraph object.
#' @param info accept two type of data: an igraph object, or KEGG.stats object obtained by \code{\link{KEGG_enrich}}.
#' @param layout_method The layout method to use: One of dot, neato, twopi, circo, and fdp. The default is dot
#' @param p.cut cut off value of BH adjust p value. Node with greater p value may not be included in the graph if it is not a parent node of any selected nodes.
#' @param plot.genes TRUE or FALSE (default). If TRUE, k genes will be shown in the graph.
#' @param edge.col edge.color
#' @param cex.label number, scaling factor for node label.
#' @param cex.node number, scaling factor for node size (shape).
#' @param cex.arrowsize number, scaling factor for arrow size.
#' @param cex.bar.label number, scaling factor for color bar indicators (p values).
#' @param bar.adj integer. Adjustment (horizontal) of color bar.
#' @param label.description TRUE/FALSE. If true, replace node label with pathway description.
#' @param ... refer to \code{\link{Rgraphviz::GraphvizAttributes}}
#' @return NULL
#' @author ZG Zhao
#' @export
plotKEGG <- function(info, layout_method='fdp', p.cut=0.05, plot.genes=FALSE, edge.col="black",
                     cex.label=1, cex.node=1, cex.arrowsize=1,
                     cex.bar.label=1, bar.adj=0,
                     label.description=TRUE, ...) {

    stopifnot(inherits(info, "KEGG.stats"))

    dtinfo <- info$stats
    gg <- kinfo2graph(rownames(dtinfo[dtinfo$p.adj < p.cut, ]), include.genes = plot.genes)
    nn <- numNodes(gg)
    nodeNames <- nodes(gg)
    dtinfo <- dtinfo[nodeNames, ]
    ## node size ~ percent.gene
    nodeSize <- dtinfo[nodeNames, "percent.test"] * 1000 * cex.node
    nodeSize <- round(sqrt(nodeSize), 0)
    ## node text
    ss1 <- nodeNames %in% keggs$ko[keggs$level == "A"]
    ss3 <- nodeNames %in% keggs$ko[keggs$level == "C"]
    ss4 <- nodeNames %in% keggs$ko[keggs$level == "D"]
    textColor <- rep("blue", nn)
    textColor[ss4] <- "transparent"
    ## labels
    if(label.description) {
        pathinfo <- keggs[!duplicated(keggs$ko), ]
        rownames(pathinfo) <- pathinfo$ko
        nodeLabel <- pathinfo[nodeNames, "description"]
    } else nodeLabel <- nodeNames

    ## label cex: more nodes, larger cex
    cex.label <- cex.label * nn^0.2 * 4
    cexLabel <- rep(cex.label, nn)
    ## NO EFFECTS !!
    ## cexLabel[ss1] <- cex.label* 3
    ## cexLabel[ss3] <- cex.label* 0.6
    ## cexLabel[ss4] <- cex.label* 0.3
    
    ## color settings for p-value
    colorPanel <- colorRampPalette(c("yellow", "orangered"))(71 - 12)
    colorPanel <- c(rep("white", 12), colorPanel)
    pvalue <- dtinfo[nodeNames, "p.adj"]
    pvalue[is.na(pvalue)] <- 1
    pvalue[pvalue < 1e-7] <- 1e-7
    pvalue <- abs(round(log10(pvalue) * 10, 0))
    nodeColor <- colorPanel[pvalue + 1]
    names(textColor) <- names(nodeColor) <- nodeNames
    names(nodeLabel) <- names(nodeSize) <- nodeNames

    nattrs <- list()
    gattrs <- Rgraphviz::getDefaultAttrs()
    gattrs$node$shape <- "circle"  ## circle, ellipse, plaintext
    gattrs$edge$color <- edge.col
    nattrs$label <- nodeLabel
    nattrs$width <- nodeSize
    nattrs$fillcolor <- nodeColor
    nattrs$fontcolor <- textColor
    gattrs$node$fontsize <- 12 * cex.label
    gattrs$edge$arrowsize <- 1 * cex.arrowsize

    mm <- matrix(rep(1, 144), nrow = 12)
    mm[1, 5:8 + bar.adj] <- 2
    layout(mm)
    Rgraphviz::plot(gg, layout_method, nodeAttrs=nattrs, attrs=gattrs, ...)
    opar <- par(no.readonly = TRUE)
    ## plot color bar
    par(mar=c(2, 2, 1, 2), xaxs="i", yaxs="i")
    xx <- barplot(rep(1, times=71),col=colorPanel,border=colorPanel, axes=FALSE)
    axis(side=1, at=as.vector(xx)[c(13, 20, 30, 40, 50, 60, 70)],
         labels=c(0.05, "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7"), cex.axis=cex.bar.label)
    box()
    par(opar)
    invisible(NULL)
}
