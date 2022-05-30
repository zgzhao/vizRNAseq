#' @name data
#' @title Useful data/variables
#' @description constants provided in the package
#' @details
#' All variables are `tibbles`:
#' - ATHGENEANNO: gene annotations of Arabidopsis
#' - ATHGENEGO: mappings from gene to KEGG indentifier.
#' - GOFULL: data parsed from go.obo file downloaded from purl.obolibrary.org.
#' - KEGGFULL: data parse from ko00001.keg file downloaded from KEGG website.
#' @author ZG Zhao
#' @examples
#' library(vizRNAseq)
#' ATHGENEANNO
#' ATHGENEGO
#' GOFULL
#' KEGGFULL
NULL


#' Helper function for extracting genes of a family (Arabidopsis).
#'
#' Refer to \code{\link{ATHGENEANNO}} variable for full list of genes.
#' @title function getGeneFamily
#' @param gene.family character, name of the family such "ARF", "IAA", "SAUR"
#' @param ignore.case TRUE/FALSE, ignore case of the gene.family name
#' @return tibble with two columns: gene_id, name
#' @author ZG Zhao
#' @export
getGeneFamily <- function(gene.family, ignore.case=FALSE) {
    annos <- ATHGENEANNO
    patt <- paste0("^", gene.family, "[0-9]+$")
    if(grepl("unix", .Platform$OS, ignore.case = T)) {
        alias <- parallel::mclapply(annos$desc, FUN=function(x){
            xdesc <- strsplit(x, "|", fixed=T)[[1]]
            xdesc <- grep(patt, xdesc, value=T, ignore.case=ignore.case)
            if(length(desc) < 1) return("")
            paste(xdesc, collapse="|")
        }, mc.cores = parallel::detectCores() -1) %>% unlist
        names(alias) <- NULL
        annos$alias <- alias
    } else {
        annos$alias <- sapply(annos$desc, FUN=function(x){
            xdesc <- strsplit(x, "|", fixed=T)[[1]]
            xdesc <- grep(patt, xdesc, value=T, ignore.case=ignore.case)
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
