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

## ## ================================================
## library(vizRNAseq)
## onto_file <- '../inst/go.obo'
## ## xrl <- 'http://purl.obolibrary.org/obo/go.obo'
## ## download.file(xrl, onto_file, method = "wget")
## GOFULL <- parseOBO(onto_file) %>%
##     mutate(space=case_when(
##                grepl('biolog', space) ~ 'BP',
##                grepl('molecular', space) ~ 'MF',
##                grepl('cellular', space) ~ 'CC'),
##            space=factor(space, levels=c('BP', 'MF', 'CC')))
## GOFULL$obsolete
## GOFULL
## save(GOFULL, file="../data/GOFULL.rda")

## ## =============================================
## library(tidyverse)
## go_tair <- '../inst/gene_association.tair.gz'
## ## xrl <- 'https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/gene_association.tair.gz'
## ## download.file(xrl, go_tair, method = "wget")

## goinfo <- tibble(ll = readLines(go_tair)) %>% 
##     filter(! grepl('^!', ll)) %>%
##     filter(grepl('\tGO:[0-9]+\t', ll)) %>% 
##     filter(grepl('\tAT[0-5CM]G[0-9]{5}\t', ll)) %>%
##     pull(ll)

## aths <- sapply(goinfo, FUN=function(x){
##     ans <- strsplit(x, '\t')[[1]]
##     gid <- grep('^AT[0-5CM]G[0-9]{5}', ans, value = T)
##     kox <- grep('^GO:[0-9]+$', ans, value = T)[1]
##     if(length(gid) < 2) gid <- c(gid, gid)
##     c(gid[1], ans[3], kox, gid[2])
## })

## aths <- t(aths)
## colnames(aths) <- c('gene_id', 'name', 'go', 'desc')
## rownames(aths) <- NULL
## aths <- as_tibble(aths) %>%
##     arrange(gene_id)

## ATHGENEGO <- aths %>% distinct(gene_id, go)
## ATHGENEANNO <- aths %>%
##     select(-go) %>% distinct(gene_id, .keep_all = TRUE)

## save(ATHGENEGO, file="../data/ATHGENEGO.rda")
## save(ATHGENEANNO, file="../data/ATHGENEANNO.rda")

## ## =============================================
## setKEGGdata <- function(d.path, force.download = FALSE) {
##     f.path <- file.path(d.path, 'ko00001.keg')
##     if(! file.exists(f.path) || force.download) {
##         if(! dir.exists(d.path)) dir.create(d.path, recursive=T)
##         xurl <- paste0("http://rest.kegg.jp/get/br:ko00001")
##         try(download.file(xurl, f.path, method = "auto"), silent = TRUE)
##     }
##     if(! file.exists(f.path) || file.size(f.path) < 10)
##         stop("File does not exist or downloaded failed.")

##     ## parse text
##     KEGGFULL <- tibble(description = readLines(f.path)) %>%
##         filter(grepl('^[A-D]', description)) %>% 
##         mutate(description = sub("^([A-D]) *", "\\1 ", description)) %>% 
##         filter(grepl("^[A-D] +\\w+ +\\w+", description)) %>%
##         mutate(description = gsub('\t.+$', '', description),
##                level = sub("^([A-D]).+$", "\\1", description),
##                ko = sub("^[A-D] +(\\w+) .*$", "\\1", toupper(description)),
##                description = sub("^[A-D] +\\w+ +(.*)$", "\\1", description),
##                description = sub(" */ *", "/", description),
##                description = sub(" \\[.+$", "", description)) %>%
##         select(level, ko, description)

##     ## ABC hierarchical
##     lvs <- KEGGFULL$level
##     nn <- nrow(KEGGFULL)
##     for(ll in LETTERS[1:3]){
##         KEGGFULL[[ll]] <- NA
##         ndx <- which(lvs == ll)
##         for(i in 1:length(ndx)){
##             ff <- ndx[i] + 1
##             KEGGFULL[[ll]][ff:nn] <- KEGGFULL$ko[ndx[i]]
##         }
##     }
##     for(i in 1:3) {
##         ss <- lvs == LETTERS[i]
##         for(j in i:3) KEGGFULL[ss, LETTERS[j]] <- NA
##     }
##     KEGGFULL <- as_tibble(KEGGFULL)
##     save(KEGGFULL, file="../data/KEGGFULL.rda")
## }

## library(tidyverse)
## setKEGGdata(d.path='../inst')

