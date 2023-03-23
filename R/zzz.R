## =============================================
## library(vizRNAseq)
## onto_file <- 'inst/go.obo'
## xrl <- 'http://purl.obolibrary.org/obo/go.obo'
## if(! file.exists(onto_file))
##     download.file(xrl, onto_file, method = "wget")

## GOFULL <- parseOBO(onto_file)
## GOFULL

## save(GOFULL, file="data/GOFULL.rda")

## ## =============================================
## KEGGFULL <- getKEGGFULL(d.path='inst')
## KEGGFULL
## save(KEGGFULL, file="data/KEGGFULL.rda")

## ## =============================================
## go_tair <- 'inst/gene_association.tair.gz'
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
##     mutate(name=sub('^AT[0-9CM]G[0-9]+\\|', '', name)) %>% 
##     select(-go) %>% distinct(gene_id, .keep_all = TRUE)


## save(ATHGENEGO, file="data/ATHGENEGO.rda")
## save(ATHGENEANNO, file="data/ATHGENEANNO.rda")
