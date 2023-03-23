#' Get latest KEGGFULL if needed.
#'
#' `KEGGFULL` is already attached in `vizRNAseq` package. You may not use this function if there is no KEGG annotation updates.
#' @title File parser: download KEGG ko00001.keg file and parse to tibble
#' @param d.path character, pathway to save KEGG brite file.
#' @return tibble with same structure of KEGGFULL
#' @author ZG Zhao
#' @export
getKEGGFULL <- function(d.path) {
    f.path <- file.path(d.path, 'ko00001.keg')
    if(! dir.exists(d.path)) dir.create(d.path, recursive=TRUE)
    xurl <- paste0("http://rest.kegg.jp/get/br:ko00001")
    if(! file.exists(f.path) || file.size(f.path) < 10)
        try(download.file(xurl, f.path, method = "auto"), silent = TRUE)
    if(! file.exists(f.path) || file.size(f.path) < 10)
        stop("File does not exist or downloaded failed.")

    ## parse text
    ans <- tibble(desc = readLines(f.path) %>% trimws) %>%
        filter(grepl('^[A-D]', desc)) %>%
        mutate(desc = sub("^([A-D]) *", "\\1 ", desc)) %>%
        filter(grepl("^[A-D] +\\w+ +\\w+", desc)) %>%
        mutate(
            level = sub("^([A-D]).+$", "\\1", desc),
            ID = sub("^[A-D] +(\\w+) .*$", "\\1", toupper(desc)),
            desc = sub("^[A-D] +\\w+ +(.*)$", "\\1", desc),
            desc = sub(" */ *", "/", desc),
            desc = sub(" \\[.+$", "", desc),
            desc = pmap_chr(list(ID, desc),  function(kk, dd) {
                sub(paste0('^', kk, '[ ,;]+'), '', dd)
            })
        ) %>%
        select(ID, level, desc)

    lvs <- ans$level
    nn <- nrow(ans)
    for(ll in LETTERS[1:3]){
        ans[[ll]] <- NA
        ndx <- which(lvs == ll)
        for(i in seq_along(ndx)){
            ff <- ndx[i] + 1
            ans[[ll]][ff:nn] <- ans$ID[ndx[i]]
        }
    }
    for(i in 1:3) {
        ss <- lvs == LETTERS[i]
        for(j in i:3) ans[ss, LETTERS[j]] <- NA
    }
    ## finalize ===========================
    as_tibble(ans) %>%
        mutate(
            ancestors=pmap(list(A, B, C), function(aa, bb, cc) na.omit(c(aa, bb, cc))),
            isa=case_when(
                is.na(A) ~ NA,
                is.na(B) ~ A,
                is.na(C) ~ B,
                TRUE ~ C)) %>%
        select(ID, ancestors, isa, desc, level) %>%
        arrange(level, ID)
}

#' Get/parse gene-KO map from file
#'
#' Download (if not existed) KEGG brite file and parse to gene-KO map.
#' @title File parser: download (from KEGG website) and parse species specific KEGG pathway data
#' @param org character, for KEGG organism name
#' @param d.path character, local path having or to save KEGG brite files. Default is `KEGG/brite` in `getwd()`.
#' @param force.download logic, default FALSE. Force to download KEGG brite file or not.
#' @return tibble with columns: gene_id, gene_name, ko and desc.
#' @author ZG Zhao
#' @export
#' @examples
#' vizRNAseq::KEGG_gene_map('ath', d.path='KEGG/brite', force.download=FALSE)
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
    tibble(gene_id=gene, gene_name=name, ko=ko) %>%
        distinct(gene_id, ko, .keep_all=TRUE)
}

#' @export
KEGG_full_map <- function(gmap) {
    onto_full_map(gmap, onto.type='ko')
}

