#' Plot heatmap with pheatmap
#'
#' Heamap plot function.
#' @title Function: plotHeat
#' @param datax data.frame. Columns "gene_id" and "agi" are required together with value columns.
#' @param gids gene_id to subset.
#' @param agi agi to subset.
#' @param alias Gene name alias.
#' @param fname File path for save tiff figure.
#' @param width Figure width.
#' @param height Figure height
#' @param fontsize Base font size
#' @param fontsize_row Font size for row labels.
#' @param fontsize_col Font size for column labels.
#' @param cut.k Number to cut.
#' @param clust.row TRUE/FALSE.
#' @param ... Further params passed to pheatmap function.
#' @return NULL
#' @author ZG Zhao
#' @export
plotHeat <- function(datax, gids=NULL, agi=NULL, alias=NULL, fname=NULL, width=800, height=200, v.cut=NULL, 
                     fontsize=9, fontsize_row=12, fontsize_col=6, cut.k=20, border=NA, 
                     clust.row = FALSE, ...){
    if (! "agi" %in% colnames(datax)) datax$agi <- datax$gene_id
    
    if (! is.null(agi) ) {
        dtx <- datax[datax$agi %in% agi, ]
        if ( is.null(alias) ) alias <- dtx$agi
        else {
            names(alias) <- agi
            alias <- alias[dtx$agi]
        }
    } else {
        if (is.null(gids)) dtx <- datax
        else dtx <- datax[datax$gene_id %in% gids, ]
        alias <- dtx$agi        
    }
    
    annos <- c("gene_id", "agi", "gene_short_name", "alias")    
    rownames(dtx) <- dtx$gene_id    
    dxx <- t(dtx[, ! colnames(dtx) %in% annos])
    if (!is.null(v.cut)) {
        xcut <- v.cut + 0.6
        dxx[dxx > xcut] <- xcut
        dxx[dxx < -xcut] <- -xcut
    }

    xmin <- min(round(dxx * 10)/10)
    xmax <- max(round(dxx * 10)/10)
    cols <- NULL
    brks <- seq(xmin, xmax, by=0.1)
    if (xmin < 0) {
        mm <- length(seq(xmin, 0, by=0.1))
        cols <- c(cols, colorRampPalette(c("black", "steelblue", "white"))(mm))
    }
    if (xmax > 0) {
        nn <- length(seq(0, xmax, by=0.1))
        cols <- c(cols, colorRampPalette(c("white", "orangeRed", "darkRed"))(nn))
    }
    cols <- colorRampPalette(cols)(length(brks) - 1)
    
    if (! is.null(fname)) tiff(fname, width=width, height=height)
    pheatmap(dxx, color = cols, labels_col = alias, cutree_cols=cut.k, cluster_rows = clust.row, fontfamily="Helvetica", 
             border_color=border, breaks=brks, ...)
    if (! is.null(fname)) dev.off()
}

#' plotHeatAP
#'
#' plotHeatAP
#' @title plotHeatAP
#' @param dtx 
#' @param agi 
#' @param gids 
#' @param fname 
#' @param width 
#' @param height 
#' @param v.cut 
#' @param ... 
#' @return NULL
#' @author ZG Zhao
#' @export
plotHeatAP <- function(dtx, agi=NULL, gids=NULL, apclust=TRUE, fname=NULL, width=400, height=600, v.cut=6, ...){
    
    if (! "agi" %in% colnames(dtx)) dtx$agi <- dtx$gene_id
    if (!is.null(gids)) dtx <- dtx[dtx$gene_id %in% gids, ]
    
    rowlab <- dtx$gene_id
    if (! is.null(agi) ) {
        dtx <- dtx[dtx$agi %in% agi, ]
        rowlab <- dtx$agi
    }
    
    if (apclust) {
        if (nrow(dtx) > 2000) stop("AP data should not exceed 2000 rows!")
        
        xdist <- function(x) { negDistMat(x) }
        xclust <- function(x) {
            aggres <- aggExCluster(negDistMat(r=2), x)
            as.hclust(aggres)
        }
    }
    
    annos <- c("gene_id", "agi", "gene_short_name", "alias")    
    rownames(dtx) <- dtx$gene_id    
    dtx <- as.matrix(dtx[, ! colnames(dtx) %in% annos])
    
    ## 颜色设置
    xcut <- v.cut + 0.6
    dtx[dtx > xcut] <- xcut
    dtx[dtx < -xcut] <- -xcut
    
    xmin <- min(round(dtx * 10)/10)
    xmax <- max(round(dtx * 10)/10)
    cols <- NULL
    brks <- seq(xmin * 10, xmax * 10, by=1)
    xn <- length(brks)
    if (xmin < 0) {
        mm <- length(seq(xmin, 0, by=0.1))
        cols <- c(cols, colorRampPalette(c("black", "steelblue", "white"))(mm))
    }
    if (xmax > 0) {
        nn <- length(seq(0, xmax, by=0.1))
        cols <- c(cols, colorRampPalette(c("white", "orangeRed", "darkRed"))(nn))
    }
    cols <- colorRampPalette(cols)(xn)
    xll <- brks %in% seq(-v.cut * 10, v.cut * 10, by=20)
    
    lmat <- rbind( c(5,3,4), c(2,1,4))
    lhei <- c(1, 7)
    lwid <- c(1, 5, 0.5)
    xbar <-  function() {
        oldpar <- par(c("mar", "mgp"))
        par(mar=c(2.5, 0.5, 2, 0.1), mgp=c(1, 0.5, 0))
        xx <- barplot(rep(1, xn), col=cols, border=cols, axes=FALSE)
        axis(1, at=xx[xll], labels=brks[xll]/10, lwd=.1)
        par(oldpar)
    }
    
    if (! is.null(fname)) tiff(fname, width=width, height=height)
    if (apclust) 
        heatmap.2(dtx, distfun = xdist, hclustfun = xclust, col=cols, labRow = rowlab,
                  lmat=lmat, lhei=lhei, lwid=lwid, key=FALSE, extrafun=xbar,
                  scale="none", density.info="none", trace="none", ...)
    else
        heatmap.2(dtx, col=cols, labRow = rowlab, lmat=lmat, lhei=lhei, lwid=lwid, key=FALSE, extrafun=xbar,
                  scale="none", density.info="none", trace="none", ...)

    if (! is.null(fname)) dev.off()
}

#' plot heatmap with ggplot2
#'
#' Select a subset of genes and plot expression heatmap.
#' @title Function: plotHeat2
#' @param dtx data.frame
#' @param gids names of gene subset for plotting.
#' @param treats Sample name labels.
#' @param cut.k groups to cut (for cutree).
#' @param fname Image file name.
#' @param width figure width
#' @param height figure height
#' @param border color for geom_tile.
#' @param del columns to omit.
#' @param flip TRUE/FALSE. flip xy
#' @param gtheme additional ggplot2 theme sets.
#' @param fill.mid mid color for heatmap.
#' @param exprs Expression or log2 fold change data. First column should be gene_id.
#' @return NULL
#' @author ZG Zhao
#' @export
plotHeat2 <- function(dtx, gids, treats=NULL, cut.k=8, fname=NULL, width=1000, height=250, v.cut=0.5, 
                     border="gray", del=NULL, flip=FALSE, gtheme=NULL,
                     fill.mid="white", fill.low="steelblue", fill.high="sienna4",
                     fill.legend=TRUE){

    if (! "agi" %in% colnames(dtx)) dtx$agi <- dtx$gene_id
    dtx <- dtx[dtx$gene_id %in% gids, ]
    rownames(dtx) <- dtx$gene_id
    
    keeps <- c("gene_id", "agi", "gene_short_name")
    keeps <- keeps[keeps %in% colnames(dtx)]
    hc <- hclust(dist(dtx[, ! colnames(dtx) %in% keeps]), method="complete")
    dtx <- dtx[hc$order, ]
    dtx <- melt(dtx, varnames = keeps, variable.name = "treatment")
    dtx$gene_id <- factor(dtx$gene_id, levels = unique(dtx$gene_id))
    dtx$value <- round(dtx$value/v.cut) * v.cut
    lbs <- dtx$agi
    names(lbs) <- dtx$gene_id    
    xlab <- function(x) lbs[x]
    if(cut.k < 9) xcols <- brewer.pal(cut.k, "RdBu")
    else xcols <- colorRampPalette(c("blue", "white", "darkgreen", "yellow", "orangered"))(cut.k)
    
    xcols <- rep(xcols, length=cut.k)[cutree(hc,k=cut.k)][hc$order]
    
    if (flip) {
        p <- ggplot(dtx, aes(x=treatment, y=gene_id, fill=value)) + gtheme  + scale_y_discrete(expand = c(0, 0)) + 
        scale_x_discrete(expand = c(0, 0))
    } else {
        p <- ggplot(dtx, aes(x=gene_id, y=treatment, fill=value)) + gtheme  + scale_y_discrete(expand = c(0, 0)) + 
        scale_x_discrete(expand = c(0, 0), labels=xlab) + 
        theme(axis.text.x = element_text(angle = 90, vjust =0.5, size=10, color=xcols))
    }
    
    p <- p + geom_tile(colour = border)
    p <- p + scale_fill_gradient2(low = fill.low, mid=fill.mid, high = fill.high, midpoint=0)
    p <- p + theme(axis.ticks = element_blank(), axis.line=element_blank(), axis.title = element_blank())

    if (! fill.legend ) p <- p + scale_fill_continuous(guide=FALSE)
    if (is.null(fname)) print(p)
    else {
        tiff(fname, width = width, height = height)
        print(p)
        dev.off()
    }
}

