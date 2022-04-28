#' Parse GO term info from obo format file.
#'
#' Parse GO term info from obo format file. 
#' @title Parse obo file
#' @param f character, file path to obo file.
#' @return tibble
#' @author ZG Zhao
#' @export
parseOBO <- function(f) {
    glist <- oboSplit(f)
    tibble(go =oboID(glist),
           desc= oboName(glist),
           space = oboNameSpace(glist),
           isa = oboIsa(glist),
           cat = oboSubset(glist),
           altid = oboAltIDs(glist),
           obsolete = oboIsObsolete(glist),
           synonym = oboSynonyms(glist)
           )
}

oboSplit <- function(f) {
    oinf <- readLines(f)
    n1 <- grep('[Term]', oinf, fixed = T)[1]
    n2 <- grep('[Typedef]', oinf, fixed = T)[1] -1
    oinf <- oinf[n1:n2]
    n1 <- grep('[Term]', oinf, fixed = T)
    n2 <- grep('^ *$', oinf)
    ans <- lapply(1:length(n1), FUN = function(i){
        aa <- n1[i]
        bb <- n2[i]
        oinf[aa:bb]
    })
    ans
}
    
## =======================================================
## single and unique item
oboID <- function(ll) {
    sapply(ll, FUN=function(x){
        xx <- grep('^id:', x, value = T)
        if(length(xx) > 0) sub('^.+(GO:[0-9]+) *$', '\\1', xx[1])
        else NA
    })
}
oboName <- function(ll) {
    sapply(ll, FUN=function(x){
        xx <- grep('^name:', x, value = T)
        if(length(xx) > 0) sub('^name: *(.+) *$', '\\1', xx[1])
        else NA
    })
}

oboNameSpace <- function(ll) {
    sapply(ll, FUN=function(x){
        xx <- grep('^namespace:', x, value = T)
        if(length(xx) > 0) sub('^namespace: *(.+) *$', '\\1', xx[1])
        else NA
    })
}
oboIsObsolete <- function(ll) {
    sapply(ll, FUN=function(x){
        any(grepl('^is_obsolete: *true *', x, ignore.case=T))
    })
}

## =======================================================
## maybe multi items
oboAltIDs <- function(ll) {
    sapply(ll, FUN=function(x){
        xx <- grep('alt_id:', x, value = T)
        if(length(xx) < 1) return(NA)
        sub('^alt_id:.*(GO:[0-9]+)[^0-9]*$', '\\1', xx)
    })
}

oboIsa <- function(ll) {
    sapply(ll, FUN=function(x){
        xx <- grep('^is_a.+(GO:)', x, value = T)
        if(length(xx) < 1) return(NA)
        sub('^is_a.+(GO:[0-9]+)[^0-9]*$', '\\1', xx)
    })
}

oboSynonyms <- function(ll) {
    sapply(ll, FUN=function(x){
        xx <- grep('synonym: *".+"', x, value = T)
        if(length(xx) < 1) return(NA)
        sub('^synonym: *"(.+)".*$', '\\1', xx)
    })
}

oboSubset <- function(ll) {
    sapply(ll, FUN=function(x){
        xx <- grep('^subset:', x, value = T)
        if(length(xx) > 0) sub('^subset: *(.+) *$', '\\1', xx)
        else 'generic'
    })
}

oboConsider<- function(ll) {
    sapply(ll, FUN=function(x){
        xx <- grep('^consider:', x, value = T)
        if(length(xx) > 0) sub('^consider: *(.+) *$', '\\1', xx)
        else NA
    })
}

