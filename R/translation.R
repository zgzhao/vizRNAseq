exAAA <- function(x) {
    sss <- as.character(x)
    ans <- sapply(sss, FUN=function(xx) {
        if(! grepl('M[^\\*]{30,}\\*', xx))
            return('X')
        xxx <- strsplit(xx, "[^\\*]\\*",)[[1]]
        xxx <- grep('M', xxx, value = T)
        if(length(xxx) < 1) return('X')
        xxx <- sub('^[^M]+', '', xxx)
        ncs <- nchar(xxx)
        xxx <- xxx[ncs == max(ncs)]
        xxx[1]
    })
    ans <- AAStringSet(ans)
    sx <- grepl('^X$', ans)
    ans[! sx]
}

#' Sequence translation
#'
#' Sequence translation
#' @title Get best/longest AA sequence from fasta file (DNA strings)
#' @param fasta_file 
#' @return AAStringSet
#' @author ZG Zhao
#' @export
bestORFs <- function(fasta_file) {
    xseqs <- readDNAStringSet(fasta_file)
    rseqs <- reverseComplement(xseqs)
    proteins <- list(
        orf1 = translate(xseqs,if.fuzzy.codon = 'X'),
        orf2 = translate(DNAStringSet(xseqs, start=2),if.fuzzy.codon = 'X'),
        orf3 = translate(DNAStringSet(xseqs, start=3),if.fuzzy.codon = 'X'),
        orf4 = translate(rseqs),
        orf5 = translate(DNAStringSet(rseqs, start=2),if.fuzzy.codon = 'X'),
        orf6 = translate(DNAStringSet(rseqs, start=3),if.fuzzy.codon = 'X')
    )
    aalist <- mclapply(1:6, FUN=function(i) {
        rex <- exAAA(proteins[[i]])
        names(rex) <- paste0(names(rex), '_orf', i)
        rex
    }, mc.cores = 6)
    rex <- unlist(AAStringSetList(aalist))
    rex[ !duplicated(rex)]
}
