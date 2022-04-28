#' @export
KEGG_gene_map <- function(org, d.path='KEGG/brite', force.download = FALSE) {
    f.path <- file.path(d.path, paste0(org, '00001.txt'))
    if(! file.exists(f.path) || force.download) {
        if(! dir.exists(d.path)) dir.create(d.path, recursive=T)
        xurl <- paste0("http://rest.kegg.jp/get/br:", org, "00001")
        try(download.file(xurl, f.path, method = "auto"),
            silent = TRUE, outFile = "kegg.log")
    }
    
    if(! file.exists(f.path) || file.size(f.path) < 10) {
        warning("File does not exist or downloaded failed.")
        invisible(NULL)
    }
    linfo <- readLines(f.path) %>%
        grep('^D.+\t.+', ., value=T) %>%
        sub('^D *', '', .)
    gene <- sub("^(\\w+) .*$", "\\1", toupper(linfo))
    name <- sub("^\\w+ +(.*)\t.+$", "\\1", linfo)
    ko <- sub("^.+\t(\\w+) .+$", "\\1", linfo)
    desc <- sub("^.+\t\\w+ +(.+)$", "\\1", linfo) %>%
        sub(" *\\[.+$", '', .)
    tibble(gene_id=gene, gene_name=name, ko=ko, desc=desc)
}

#' @export
KEGG_full_map <- function(gene2ko){
    kxx <- c('09100', '09120', '09130', '09140', '09150')
    KEGGFULL %>%
        filter(ko %in% {{kxx}} | A %in% {{kxx}}) %>% 
        left_join(
               gene2ko %>% group_by(ko) %>%
               summarise(gene = list(unique(gene_id))),
               by = 'ko') %>%
        select(ko, A, B, C, gene) %>% 
        unnest(gene) %>%
        pivot_longer(cols=-gene, values_to = 'ko') %>%
        select(ko, gene) %>% distinct
}
