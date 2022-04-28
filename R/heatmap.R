
#' @export
vheatmap <- function(mat, vcut=4, border=TRUE, ...) {
    vcut <- abs(vcut)
    brks <- seq(-vcut, vcut, by=0.1)
    nn <- (length(brks) -1)/2
    col1 <- colorRampPalette(c("black", "steelblue", "white"))(nn)
    col2 <- colorRampPalette(c("white", "orangeRed", "darkRed"))(nn)
    cols <- c(col1, col2)
    bcol <- if(border) 'gray60' else 'transparent'
    pheatmap(mat, color=cols, breaks = brks, border_color = bcol, ...)
}

#' Plot heatmap with pheatmap
#'
#' Heamap plot function.
#' @title Function: plotHeat
#' @param datax data.frame. Columns "gene_id" and "agi" are required together with value columns.
#' @param gids gene_id to subset.
#' @param agi agi to subset.
#' @param alias Gene name alias.
#' @param trans.axis TRUE/FALSE (default)
#' @param v.cut cut value
#' @param cut.k Number to cut.
#' @param border grid border color
#' @param clust.row TRUE/FALSE.
#' @param clust.col TRUE/FALSE.
#' @param ... Further params passed to pheatmap function.
#' @return NULL
#' @author ZG Zhao
#' @export
plotHeat <- function(datax, gids=NULL, agi=NULL, alias=NULL,
                     trans.axis=FALSE, v.cut=4, cut.k=20, border=NA,
                     clust.row = FALSE, clust.col=TRUE, ...){
    if (! "agi" %in% colnames(datax)) datax$agi <- datax$gene_id

    ## --------------------
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
        if(is.null(alias)) alias <- dtx$agi
    }

    ## --------------------
    annos <- c("gene_id", "agi", "gene_short_name", "alias")    
    rownames(dtx) <- dtx$gene_id    
    dxx <- t(dtx[, ! colnames(dtx) %in% annos])
    if (!is.null(v.cut)) {
        xcut <- v.cut + 0.6
        dxx[dxx > xcut] <- xcut
        dxx[dxx < -xcut] <- -xcut
    } else xcut <- ceiling(max(dxx))

    ## --------------------
    xmin <- min(round(dxx * 10)/10)
    xmax <- max(round(dxx * 10)/10)
    col1 <- col2 <- NULL
    xrange <- brks <- seq(xmin, xmax, by=0.1)
    if(xmin > -xcut) xrange <- c(seq(-xcut, xmin - 0.1, by=0.1), xrange)
    if(xmax < xcut) xrange <- c(xrange, seq(xmax + 0.1, xcut, by=0.1))
    nn <- floor(length(xrange) / 2)
    col1 <- colorRampPalette(c("black", "steelblue", "white"))(nn)
    col2 <- colorRampPalette(c("white", "orangeRed", "darkRed"))(nn + 1)
    cols <- c(col1, col2)
    cols <- cols[xrange %in% brks]
    cols <- colorRampPalette(cols)(length(brks))
    ## --------------------
    if(trans.axis) {
        dxx <- t(dxx)
        pheatmap(dxx, color = cols, labels_row = alias, cutree_rows=cut.k,
                 cluster_rows = clust.row, cluster_cols=clust.col, fontfamily="Helvetica",
                 border_color=border, breaks=brks, ...)
    } else {
        pheatmap(dxx, color = cols, labels_col = alias, cutree_cols=cut.k,
                 cluster_rows = clust.row, cluster_cols=clust.col, fontfamily="Helvetica",
                 border_color=border, breaks=brks, ...)
    }
}

#' plot heatmap with ggplot2
#'
#' Select a subset of genes and plot expression heatmap.
#' @title Function: plotGHeat
#' @param dtx data.frame
#' @param gids names of gene subset for plotting.
#' @param treats Sample name labels.
#' @param cut.k groups to cut (for cutree).
#' @param border color for geom_tile.
#' @param del columns to omit.
#' @param flip TRUE/FALSE. flip xy
#' @param gtheme additional ggplot2 theme sets.
#' @param fill.mid mid color for heatmap.
#' @param exprs Expression or log2 fold change data. First column should be gene_id.
#' @return NULL
#' @author ZG Zhao
#' @export
plotGHeat <- function(dtx, gids, treats=NULL, cut.k=8, v.cut=0.5,
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
    print(p)
}

.trans_exprs <- function(x) {
    x[x < 1] <- 1
    log10(x)
}

#' @export
exprsHeatmap <- function(exprs_df,
                          d_cols=NULL, lb_col='gene_id',
                          g_ids=NULL, trans_axis=FALSE,
                          ...){
    if(! 'gene_id' %in% colnames(exprs_df))
        stop('Must contain `gene_id` column!')
    if(any(duplicated(exprs_df)))
        stop('Duplicated gene id is not allowed!')

    if(! is.null(g_ids))
        exprs_df <- filter(exprs_df, gene_id %in% g_ids)
    if(is.null(d_cols)) exprs_mx <- select(exprs_df, where(is.numeric))
    else exprs_mx <- exprs_df[, d_cols]
    
    exprs_mx <- as.matrix(exprs_mx) %>% .trans_exprs
    ss <- apply(exprs_mx, 1, FUN=function(x) !all(x==0))
    exprs_mx <- exprs_mx[ss, ]
    exprs_df <- exprs_df[ss, ]
    if(trans_axis) {
        exprs_mx <- t(exprs_mx)
        pheatmap(exprs_mx, scale='row',
                 cluster_rows = FALSE, cluster_cols = TRUE,
                 labels_col=exprs_df[[lb_col]], ...)
    } else {
        pheatmap(exprs_mx, scale='row',
                 cluster_rows = TRUE, cluster_cols = FALSE,
                 labels_row=exprs_df[[lb_col]], angle_col=0, ...)
    }
}
