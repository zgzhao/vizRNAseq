
.setDEGtype <- function(v, p, vc, pc) {
    v[is.na(v)] <- 0
    p[is.na(p)] <- 0
    ans <- rep(0, length(v))
    ans[abs(v) > vc & p < pc] <- -1
    ans[v > vc & p < pc] <- 1
    ans
}

#' @export
DEG_detect <- function(dfv, dfp, v.cut, p.cut) {
    dfv <- arrange(dfv, gene_id)
    dfp <- arrange(dfp, gene_id)
    gene1 <- dfv$gene_id
    gene2 <- dfp$gene_id
    if(! all(gene1 == gene2) || !all(names(dfv) == names(dfp)))
        stop('LogFC and Pvalue tables not matched!')
    for(i in 2:ncol(dfv)) {
        dfv[[i]] <- .setDEGtype(dfv[[i]], dfp[[i]], v.cut, p.cut)
    }
    dfv
}

#' @export
DEG_filter <- function(dfv, dfp, v.cut, p.cut,
                       f.mode=c("any", "all"),
                       compare=c("all", "up", "down")) {
    f.mode <- f.mode[1]
    compare <- compare[1]
    dfs <- DEG_detect(dfv, dfp, v.cut, p.cut)
    ss <- apply(dfs[, -1], 1, FUN=function(x) all(x == 0))
    dfs <- filter(dfs, !ss)
    genes <- dfs$gene_id
    
    xfun <- if(f.mode == "any") any else all
    if(compare == "up") {
        ss <- apply(dfs[, -1], 1, FUN=function(x) xfun(x > 0))
        genes <- genes[ss]
    } else if(compare == "down") {
        ss <- apply(dfs[, -1], 1, FUN=function(x) xfun(x < 0))
        genes <- genes[ss]
    }
    dfv %>% filter(gene_id %in% {{genes}})
}
