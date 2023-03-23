#' TODO
#'
#' TODO
#' @aliases KEGG_full_map GO_full_map
#' @title convert partial gene annotation map to full KEGG/GO tree map
#' @param gmap data.frame with two columns: gene_id and ko (or go)
#' @param onto.type character scalar, accept 'ko' or 'go' only
#' @return tibble with two columns: ID and genes
#' @author ZG Zhao
#' @export
onto_full_map <- function(gmap, onto.type=c('ko', 'go')) {
    if (onto.type[1] == 'ko') {
        onto_full <-  KEGGFULL
        scol <- 'ko'
    } else {
        onto_full <-  GOFULL
        scol <- 'go'
    }
    missing_cols <- setdiff(c('gene_id', scol), colnames(gmap))
    if (length(missing_cols) > 0) {
        msn <- paste(missing_cols, collapse = ', ')
        stop(paste('Missing columns in gmap:', msn))
    }
    gmap <- gmap[, c(scol, 'gene_id')]
    colnames(gmap) <- c('ID', 'genes')
    ## formate alternative GO ids
    if (scol == 'go') {
        gmap <- gmap %>%
            left_join(
                GOFULL %>%
                select(ID, altid) %>%
                rename(xid=ID) %>%
                rename(ID=altid) %>% 
                unnest(ID),
                by = 'ID') %>%
            mutate(ID=ifelse(is.na(xid), ID, xid)) %>%
            select(-xid)
    }
    
    gmap %>% left_join(onto_full, by = 'ID') %>%
        select(ancestors, genes) %>%
        unnest(ancestors) %>%
        rename(ID=ancestors) %>%
        bind_rows(gmap) %>% 
        group_by(ID) %>%
        summarise(genes=list(unique(unlist(genes)))) %>%
        mutate(N=map_int(genes, length)) %>%
        filter(N > 0) %>%
        select(-N) %>%
        arrange(ID)
}


#' Helper function for extracting genes of a family (Arabidopsis).
#'
#' Refer to \code{\link{ATHGENEANNO}} variable for full list of genes.
#' @title Arabidopsi: get genes by family keyword
#' @param keyword character, name of the family such "ARF", "IAA", "SAUR"
#' @param ignore.case TRUE/FALSE, ignore case of the keyword
#' @return tibble with two columns: gene_id, name
#' @author ZG Zhao
#' @export
ath_gene_family <- function(keyword, ignore.case=FALSE) {
    annos <- ATHGENEANNO
    patt <- paste0("^", keyword, "[0-9]+$")
    if(grepl("unix", .Platform$OS, ignore.case = T)) {
        alias <- parallel::mclapply(annos$desc, FUN=function(x) {
            xdesc <- strsplit(x, "|", fixed=TRUE)[[1]]
            xdesc <- grep(patt, xdesc, value=TRUE, ignore.case=ignore.case)
            if(length(desc) < 1) return("")
            paste(xdesc, collapse="|")
        }, mc.cores = parallel::detectCores() -1) %>% unlist
        names(alias) <- NULL
        annos$alias <- alias
    } else {
        annos$alias <- sapply(annos$desc, FUN=function(x) {
            xdesc <- strsplit(x, "|", fixed=TRUE)[[1]]
            xdesc <- grep(patt, xdesc, value=TRUE, ignore.case=ignore.case)
            if(length(desc) < 1) return("")
            paste(xdesc, collapse="|")
        })
    }
    annos <- annos[grepl(patt, annos$name) |
                   annos$alias != "", ]
    ss <- ! grepl(patt, annos$name)
    annos$name[ss] <- annos$alias[ss]
    annos %>% select(gene_id, name) %>%
        mutate(ndx = sub("^[^0-9]+([0-9]+).*$", "\\1", name)) %>%
        arrange(as.integer(ndx)) %>%
        select(-ndx)
}


