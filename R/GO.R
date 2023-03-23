#' @export
GO_full_map <- function(gmap) {
    onto_full_map(gmap, onto.type='go')
}

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
    ans <- tibble(ID=oboID(glist),
                  isa=oboIsa(glist),
           desc=oboName(glist),
           space=oboNameSpace(glist),
           cat=oboSubset(glist),
           altid=oboAltIDs(glist),
           obsolete=oboIsObsolete(glist),
           synonym=oboSynonyms(glist)) %>%
        mutate(
            isa=map(isa, function(x) na.omit(as.character(x))),
            altid=map(altid, function(x) na.omit(as.character(x))),
            synonym=map(synonym, function(x) na.omit(as.character(x))),
            space=case_when(
                grepl('biolog', space) ~ 'BP',
                grepl('molecular', space) ~ 'MF',
                grepl('cellular', space) ~ 'CC'))
    
    res <- NULL
    for (cc in c('CC', 'MF', 'BP')) {
        anx <- filter(ans, space=={{cc}})
        xdist <- anx %>% 
            select(ID, isa) %>%
            unnest(isa) %>% 
            select(isa, ID) %>% 
            igraph::graph.data.frame() %>% 
            igraph::distances(mode = 'in')
        res <- as_tibble(xdist) %>%
            mutate(ID=rownames(xdist)) %>%
            pivot_longer(-ID) %>% 
            filter(value > 0) %>%
            filter(! is.infinite(value)) %>%
            group_by(ID) %>%
            summarise(connected=list(c(name))) %>%
            right_join(anx, by='ID') %>%
            bind_rows(res)
    }
    res
}

oboSplit <- function(f) {
    oinf <- readLines(f)
    n1 <- grep('[Term]', oinf, fixed = T)[1]
    n2 <- grep('[Typedef]', oinf, fixed = T)[1] -1
    oinf <- oinf[n1:n2]
    n1 <- grep('[Term]', oinf, fixed = T)
    n2 <- grep('^ *$', oinf)
    ans <- lapply(1:length(n1), FUN = function(i) {
        aa <- n1[i]
        bb <- n2[i]
        oinf[aa:bb]
    })
    ans
}
    
## =======================================================
## single and unique item
oboID <- function(ll) {
    sapply(ll, FUN=function(x) {
        xx <- grep('^id:', x, value = T)
        if(length(xx) > 0) sub('^.+(GO:[0-9]+) *$', '\\1', xx[1])
        else NA
    })
}
oboName <- function(ll) {
    sapply(ll, FUN=function(x) {
        xx <- grep('^name:', x, value = T)
        if(length(xx) > 0) sub('^name: *(.+) *$', '\\1', xx[1])
        else NA
    })
}

oboNameSpace <- function(ll) {
    sapply(ll, FUN=function(x) {
        xx <- grep('^namespace:', x, value = T)
        if(length(xx) > 0) sub('^namespace: *(.+) *$', '\\1', xx[1])
        else NA
    })
}
oboIsObsolete <- function(ll) {
    sapply(ll, FUN=function(x) {
        any(grepl('^is_obsolete: *true *', x, ignore.case=T))
    })
}

## =======================================================
## maybe multi items
oboAltIDs <- function(ll) {
    sapply(ll, FUN=function(x) {
        xx <- grep('alt_id:', x, value = TRUE)
        if(length(xx) < 1) return(NA)
        xx <- strsplit(xx, ' +')[[1]]
        grep('GO:[0-9]+', xx, value=TRUE)
    })
}

oboIsa <- function(ll) {
    lapply(ll, FUN=function(x) {
        xx <- grep('^is_a.+(GO:)', x, value = TRUE)
        if(length(xx) < 1) return(NA)
        xx <- strsplit(xx, ' +')[[1]]
        grep('GO:[0-9]+', xx, value=TRUE)
    })
}

oboSynonyms <- function(ll) {
    sapply(ll, FUN=function(x) {
        xx <- grep('synonym: *".+"', x, value = TRUE)
        if(length(xx) < 1) return(NA)
        xx <- strsplit(xx, ' +')[[1]]
        grep('GO:[0-9]+', xx, value=TRUE)
    })
}

oboSubset <- function(ll) {
    sapply(ll, FUN=function(x) {
        xx <- grep('^subset:', x, value = TRUE)
        if(length(xx) > 0) sub('^subset: *(.+) *$', '\\1', xx)
        else 'generic'
    })
}

oboConsider<- function(ll) {
    sapply(ll, FUN=function(x) {
        xx <- grep('^consider:', x, value = TRUE)
        if(length(xx) > 0) sub('^consider: *(.+) *$', '\\1', xx)
        else NA
    })
}
