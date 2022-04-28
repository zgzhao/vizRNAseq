.venn3 <- function(dtx, fname, cnames=NULL, fig.size=600, adj.cat=0, adj.pos=0, adj.mar=0, cex.num=0, cex.cat=0) {
    col3 <- c("red", "blue", "green")
    if (is.null(cnames)) cnames <- names(dtx)
    venn.diagram(x = dtx, height = fig.size, width = fig.size,
                 filename=fname,
                 imagetype="tiff", 
                 cex=0.3 + cex.num, cat.cex=0.3 + cex.cat, 
                 col="gray50", lwd=0.3, fill=col3, alpha=rep(0.3, 3),
                 category.names = cnames,
                 cat.pos=c(330 - adj.pos, 30 + adj.pos, 180),
                 cat.dist=c(0.1, 0.1, 0.05) + adj.cat,
                 cat.col=col3,
                 margin=0.1 + adj.mar
                 )
}

.venn4 <- function(dtx, fname, cnames=NULL, fig.size=600, adj.cat=0, adj.pos=0, adj.mar=0, cex.num=0, cex.cat=0) {
    col4 <- c("#483D8B", "#DB7093", "#FFD700", "#006400")
    if (is.null(cnames)) cnames <- names(dtx)
    venn.diagram(x = dtx, height = fig.size, width = fig.size,
                 filename=fname,
                 imagetype="tiff", 
                 cex=0.3 + cex.num, cat.cex=0.3 + cex.cat, 
                 col="gray50", lwd=0.3, fill=col4, alpha=rep(0.5, 4),
                 category.names = cnames,
                 cat.pos=c(180 - adj.pos, 180 + adj.pos, 310 - adj.pos, 50 + adj.pos),
                 cat.dist=c(0.4, 0.4, 0.18, 0.18) + adj.cat,
                 cat.col=col4,
                 margin=0.05 + adj.mar
                 )
}

#' plot venn diagram.
#'
#' Plot venn diagram with list data of length 3 or 4.
#' @title plotVenn function
#' @param dtx list of length 3 or 4.
#' @param fname String, file name to save.
#' @param cnames String vector for category names.
#' @param fig.size Numeric. Figure size.
#' @param adj.cat 
#' @param adj.pos 
#' @param adj.mar 
#' @param cex.num 
#' @param cex.cat 
#' @return NULL
#' @author ZG Zhao
#' @export
plotVenn <- function(dtx, fname, cnames=NULL, fig.size=600, adj.cat=0, adj.pos=0, adj.mar=0, cex.num=0, cex.cat=0) {
    nn <- length(dtx)
    if ( nn == 3)
        .venn3(dtx, fname, cnames=cnames, fig.size=fig.size, adj.cat=adj.cat, adj.pos=adj.pos, adj.mar=adj.mar,
               cex.num=cex.num, cex.cat=cex.cat)
    else if ( nn == 4)
        .venn4(dtx, fname, cnames=cnames, fig.size=fig.size, adj.cat=adj.cat, adj.pos=adj.pos, adj.mar=adj.mar,
               cex.num=cex.num, cex.cat=cex.cat)
    else stop("Input data must be a list of length 3 or 4.")
    logs <- list.files(dirname(fname), "\\.log$", full=T)
    if(length(logs) > 0) file.remove(logs)
    return(invisible(NULL))
}
