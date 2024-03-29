#' theme_sci
#'
#' ggplot2 theme_sci, a clear theme for figures of scientific papers.
#' @title theme_sci ggplot2 theme.
#' @param ... args passed to theme_bw().
#' @return ggplot2 theme set
#' @author ZG Zhao
#' @export
theme_sci <- function(base_size = 12, base_family = "Helvetica", ...) {
    theme_bw(base_size, base_family, ...) +
    theme(panel.grid=element_blank(),    #去掉网格线
          legend.title=element_blank(),  #去掉图例title
          legend.key=element_rect(color="white",fill="white")
      )
}

#' plotFPKM barchart
#'
#' plotFPKM barchart
#' @title plotFPKM function.
#' @param dtx data.frame
#' @param agi agi for data subset
#' @param alias gene name alias
#' @param gids gene_id filters.
#' @param filename image filename
#' @param width figure width
#' @param height figure height
#' @param sig.size size of significant labels
#' @param y.expand expand factor of y upper limit
#' @param nrow 
#' @param ncol 
#' @param cols fill colors of bars
#' @param gtheme additional theme set for ggplot object.
#' @return NULL
#' @author ZG Zhao
#' @export
plotFPKM <- function(dtx, agi=NULL, alias=NULL, gids=NULL,
                     filename=NULL, width=400, height=300, sig.size=12, y.expand=1.05, nrow=NULL, ncol=NULL, 
                     cols=c("gray80", "gray30"), gtheme=NULL){

    if (! "agi" %in% colnames(dtx)) dtx$agi <- dtx$gene_id
    if (!is.null(gids)) dtx <- dtx[dtx$gene_id %in% gids, ]
    
    if ( !is.null(agi) ) {
        dtx <- dtx[dtx$agi %in% agi, ]
        dtx$agi <- factor(dtx$agi, levels=agi)
        if ( !is.null(alias) ) {
            names(alias) <- agi
            dtx$alias <- alias[dtx$agi]            
        }
        dtx <- dtx[order(as.integer(dtx$agi)), ]
    } else dtx$alias <- dtx$agi
    alias <- dtx$alias
    names(alias) <- dtx$gene_id
    dtx$gene_id <- factor(dtx$gene_id, levels = unique(dtx$gene_id))
    
    
    if (nrow(dtx) < 1) stop("No data found for your agi.")
    
    xlabeller <- function(labels, multi_line = FALSE){
        labels <- label_value(labels, multi_line = multi_line)
        lapply(labels, function(values) {
            values <- paste0("list(`", alias[values], "`)")
            lapply(values, function(expr) c(parse(text = expr)))
        })
    }
    
    p <- ggplot(dtx, aes(x=treatment, y=mean, fill=tr.level)) + 
    labs(x = "", y = "Expression (FPKM)") +
    geom_bar(stat="identity", position = position_dodge(), color="gray30") +
    geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd), position=position_dodge(0.9), width=.5, colour="gray30") +
    geom_text(aes(y=(mean + sd) * y.expand, label="")) + 
    geom_text(aes(y=mean + sd, label=sigs), position=position_dodge(0.9), vjust=-0.2, size=sig.size) + 
    facet_wrap(~gene_id, scales="free_y", labeller=xlabeller, nrow=nrow, ncol=ncol) + scale_fill_manual(values = cols) +
    scale_y_continuous(expand = c(0, 0))
    
    if (! is.null(gtheme)) p <- p + gtheme    
    if (is.null(filename)) {
        print(p)
    } else {
        tiff(filename, width=width, height=height)
        print(p)
        dev.off()
    }
}

#' Plot expression diagram.
#'
#' Plot expression diagram of seleted genes
#' @title Function: plotChart
#' @param exprs data.frame: rownames should be gene names and columns should be fold change data for color value.
#' @param pval data.frame: rownames should be gene names and columns should be p value.
#' @param genes character vecters for gene number subset.
#' @param alias gene name alias.
#' @param samples sample name alias, should be same length of exprs and pval.
#' @param sig.size integer. significant label size.
#' @param gtheme theme object passed to ggplot2
#' @param del column names to be omitted from whole data.
#' @param num TRUE/FALSE, whether to should p values and significant labels.
#' @param flip TRUE/FALSE, flip xy axis.
#' @param clust TRUE/FALSE, cluster and reorder genes.
#' @param fname Character: string for image (tiff) file name.
#' @param width figure width
#' @param height figure height
#' @param ... not implemented currently.
#' @return NULL
#' @author ZG Zhao
#' @export
plotChart <- function(dtx, dts, agi=NULL, alias=NULL, samples=NULL, id.label=FALSE, gtheme=NULL, sig.size=5, clust=FALSE, 
                      del=NULL, num=FALSE, flip=FALSE, fname=NULL, width=400, height=300, ...){

    if (! "agi" %in% colnames(dtx)) dtx$agi <- dtx$gene_id
    if ( is.null(agi) )  stop("agi must be set.")
    
    keeps <- c("gene_id", "agi", "gene_short_name", "alias")
    
    dtx <- dtx[dtx$agi %in% agi, ]
    dts <- dts[dts$agi %in% agi, ]
    if ( !is.null(alias) ) {
        names(alias) <- agi
        dtx$alias <- alias[dtx$agi]            
    } else dtx$alias <- dtx$agi

    rownames(dtx) <- dtx$gene_id
    rownames(dts) <- dts$gene_id
    
    if (clust) {
        lfc2 <- dtx[, ! colnames(dtx) %in% keeps]
        hc <- hclust(dist(lfc2), method="complete")
        dtx <- dtx[hc$order, ]
        dts <- dts[hc$order, ]
    } else {
        dtx$agi <- factor(dtx$agi, levels=agi)
        ndx <- order(dtx$agi)
        dtx <- dtx[ndx, ]
    }

    keeps <- colnames(dtx)[colnames(dtx) %in% keeps]    
    dtx <- melt(dtx, varnames = keeps, variable.name = "lbs", value.name = "lfc")
    dts <- melt(dts, varnames = keeps, variable.name = "lbs", value.name = "pval")
    dtx$num <- round(dtx$lfc, 2)
    dtx$num[dts$pval < 0.05] <- paste(dtx$num[dts$pval < 0.05], "*")
    dtx$num[dts$pval < 0.01] <- paste0(dtx$num[dts$pval < 0.01], "*")

    dtx$gene_id <- factor(dtx$gene_id, levels=unique(dtx$gene_id))
    aliax <- dtx$alias
    names(aliax) <- dtx$gene_id    
    xlab <- if (id.label) function(x) x else function(x) aliax[x]
    if (!is.null(del)) dtx <- dtx[! dtx$gene_id %in% del, ]
    
    if (flip) {
        p <- ggplot(dtx, aes(x=lbs, y=gene_id, fill=lfc)) + gtheme + 
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_discrete(labels=xlab, expand = c(0, 0))
    } else {
        p <- ggplot(dtx, aes(x=gene_id, y=lbs, fill=lfc)) + gtheme + 
        scale_x_discrete(labels=xlab, expand = c(0, 0)) + 
        scale_y_discrete(expand = c(0, 0))
    }
    
    p <- p  + geom_tile(colour = "gray")
    if (num) p <- p + geom_text(aes(label=num), size=sig.size, colour="black")
    
    p <- p + labs(x = "", y = "") +
    scale_fill_gradient2(low = "steelblue", mid="white", high = "orangeRed", name=NULL, midpoint=0)
    
    tiff(fname, width=width, height=height)
    print(p)
    dev.off()
}

