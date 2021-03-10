#' readGIDs
#'
#' readGIDs read in gene_id or agi
#' @title Function readGIDs
#' @param files expression file
#' @param reg regulation type: all, up or down
#' @param agi TRUE/FALSE. Return agi or gene_id instead.
#' @param fc fold change filter
#' @param labels String vector for data labels.
#' @param FPKM.cut Numeric. FPKM cut value.
#' @param q.cut Numeric. q_value cut threshold.
#' @param ensemble String vector.
#' @param cc column to return
#' @author ZG Zhao
#' @export
readGIDs <- function(files, reg=c("all", "up", "down"), agi=FALSE, fc=1.5,
                     labels=NULL, FPKM.cut=0, q.cut=1, p.cut=1, ensemble=NULL) {
    reg <- reg[1]
    xcc <- if (agi) "agi" else "gene_id"

    results <- list()
    for (ff in files) {
        dtx <- read.csv(ff)
        dtx <- dtx[dtx$max_FPKM > FPKM.cut & dtx$q_value < q.cut & dtx$p_value < p.cut, ]
        
        ct1 <- dtx$log2FC_adj > log2(fc)
        ct2 <-  dtx$log2FC_adj < -log2(fc)

        if (reg == "up") xsel <-  ct1
        else if (reg == "down") xsel <- ct2
        else xsel <- ct1 | ct2
        
        dtx <- dtx[xsel, xcc]
        dtx <- unique(dtx[!is.na(dtx)])
        
        if (! is.null(ensemble) ) dtx <- intersect(ensemble, dtx)
        results <- c(results, list(dtx))
    }

    n <- length(files)
    if ( n < 2) results <- results[[1]]
    else {
        if (is.null(labels)) names(results) <- LETTERS[1:n]
        else names(results) <- labels
    }
    results
}

#' Read in and merge data.
#'
#' Read data from cuffdiff result file and merge.
#' @title mergeRead function
#' @param files filenames
#' @param d.col data column to merge, one column only.
#' @param treats Character vector of same length of filenames for names of treatments.
#' @param sels NULL. Character vector of gene_id to select.
#' @param keep TRUE/FALSE. keep the merge.by column in result data.frame.
#' @param ... args passed to read.csv
#' @return data.frame
#' @author ZG Zhao
#' @export
mergeRead <- function(files, d.col="log2FC_adj", treats=NULL, sels=NULL, keep=TRUE, ...) {
    xcol <- c("gene_id", "agi", "gene_short_name")
    datax <- read.csv(files[1], ...)
    xcol <- xcol[xcol %in% colnames(datax)]
    datax <- datax[, colnames(datax) %in% c(xcol, d.col)]
    
    if ( !is.null(sels) ) datax[datax$gene_id %in% sels, ]
    
    for ( ff in files[-1]) {
        dtx <- read.csv(ff, ...)[, d.col]
        if ( !is.null(sels) ) dtx <- dtx[dtx$gene_id %in% sels, ]
        datax <- cbind(datax, dtx)
    }
    
    if (length(treats) == ncol(datax) - length(xcol)) colnames(datax) <- c(xcol, treats)
    
    datax
}

meltData <- function(dt1, dt2, value.name="value", treats=NULL, tr.levels=NULL) {
    if ( ! all(nrow(dt1) == nrow(dt2)) ) stop("Row length of dataframes must be identical.")
    keeps <- c("gene_id", "agi", "gene_short_name")
    keeps <- keeps[keeps %in% colnames(dt1)]
    
    if (is.null(tr.levels)) tr.levels <- c("CK", "T")
    if (is.null(treats)) treats <- rev(rev(colnames(dt1))[-1])

    dtx <- rbind(dt1, dt2)
    dtx$tr.level <- rep(tr.levels, each=nrow(dt1))
    
    dtx <- melt(dtx, varnames = c(keeps, "tr.level"), variable.name = "treatment", value.name = value.name)
    dtx$treatment <- factor(dtx$treatment, levels = treats)
    dtx$tr.level <- factor(dtx$tr.level, levels = tr.levels)
    colnames(dtx) <- c(keeps, "tr.level", "treatment", value.name)
    dtx
 }

#' Read FPKM data
#'
#' Read FPKM data from files, with or without stdev and significant tr.levels.
#' @title Function readFPKM
#' @param xfiles filename
#' @param treats Character vector of same length of filenames for names of treatments.
#' @param tr.levels Character vector of 2 for treatment levels: control and treatment.
#' @param std TRUE/FALSE. combine stdev or not
#' @param sig TRUE/FALSE. combine significance label or not.
#' @param xpkm return normalized FPKM values
#' @param control If all samples share a common control, provide its name.
#' @return data.frame
#' @author ZG Zhao
#' @export
readFPKM <- function(xfiles, treats, tr.levels=NULL, std=TRUE, sig=TRUE, xpkm=FALSE, control=""){
    if (xpkm) {
        vv1 <- "XPKM1"; vv2 <- "XPKM2"; ss1 <- "xsd1"; ss2 <- "xsd2"
    } else {
        vv1 <- "value_1"; vv2 <- "value_2"; ss1 <- "stdev1"; ss2 <- "stdev2"
    }
    
    fpkm1 <- mergeRead(xfiles, d.col=vv1, treats = treats)
    fpkm2 <- mergeRead(xfiles, d.col=vv2, treats = treats)
    dtx <- meltData(fpkm1, fpkm2, "mean", treats, tr.levels)

    if (std) {
        stdx1 <- mergeRead(xfiles, d.col=ss1, treats = treats)
        stdx2 <- mergeRead(xfiles, d.col=ss2, treats = treats)
        dtx$sd <- meltData(stdx1, stdx2, "error", treats, tr.levels)$error
    }
    
    if (sig) {
        datap <- mergeRead(xfiles, d.col="p_value", treats = treats)
        pv <- meltData(datap, datap, "pval", treats, tr.levels)$pval
        sigs <- rep("", length(pv))
        sigs[pv < 0.05] <- "*"
        sigs[pv < 0.01] <- "**"
        dtx$sigs <- sigs
        ck <- if(is.null(tr.levels)) "CK" else tr.levels[1]
        dtx$sigs[dtx$tr.level == ck] <- ""
    }

    if ( control != "") {
        lv1 <- levels(dtx$tr.level)[1]
        dtx$treatment <- as.character(dtx$treatment)
        dxx <- dtx[dtx$tr.level==lv1 & dtx$treatment==treats[1], ]
        dxx$treatment <- control
        dxx <- rbind(dxx, dtx[dtx$tr.level != lv1, ])
        dxx$treatment <- factor(dxx$treatment, levels = c(control, treats))
        dtx <- dxx
    }
    
    dtx
}
