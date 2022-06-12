mseq_len <- function(ss, ee) {
    xx <- IRanges(start=ss, end = ee)
    width(IRanges::reduce(xx))
}

#' parse blast result
#'
#' parse blast result.
#' @title Get best matches from blast results.
#' @param blast_result BLAST out put file in fmt6 format
#' @param blast_db BLAST db file (fasta) when you do blast 
#' @param prot TRUE/FALSE. Protein blast db (not query sequences).
#' @param force.perl TRUE/FALSE, force to use perl program for sequence processing.
#' @return tibble
#' @author ZG Zhao
#' @export
bestMatches <- function(blast_result, blast_db, prot=FALSE, force.perl=TRUE) {
    cat("Running, please wait a few minutes ...\n")
    dtx <- read.table(blast_result) %>% as_tibble
    colnames(dtx) <- c('gene_id', 'target', 'p_indentical',
                      'n_align', 'n_mismatch', 'gap_open',
                      'qstart', 'qend',
                      'sstart', 'send',
                      'evalue', 'bit_score')
    ss <- apply(dtx[, 9:10], 1, min)
    ee <- apply(dtx[, 9:10], 1, max)
    ans <- dtx %>%
        mutate(sstart = {{ss}}, send={{ee}}) %>% 
        group_by(gene_id, target) %>%
        summarise(n_match = mseq_len(sstart, send),
                  s_mean = mean(bit_score),
                  e_mean = mean(evalue)) %>%
        left_join(dtx, by=c("gene_id", "target"))

    xperl <- file.size(blast_db) > 100000000 ## 100M
    if(xperl || force.perl) {
        pperl <- file.path(path.package("vizRNAseq"), "Perl", "filter.seqs.pl")
        xblast <- tempfile()
        xfasta <- tempfile()
        write.table(ans, xblast, row.names=F, col.names=F, quote=F, sep="\t")
        cmd <- paste("perl", pperl, xblast, blast_db, xfasta)
        system(cmd)
        if(! file.exists(xfasta))
            stop("Perl program failed!\n")
        
        ss <- if(prot) readAAStringSet(xfasta) else readDNAStringSet(xfasta)
        on.exit({
            if(file.exists(xblast)) file.remove(xblast)
            if(file.exists(xfasta)) file.remove(xfasta)
        })
    } else {
        ss <- if(prot) readAAStringSet(blast_db) else readDNAStringSet(blast_db)
    }
    xnames <- sub(' .+$', '', names(ss))
    tinfo <- tibble(target=xnames,
                    n_target=Biostrings::nchar(ss))

    ans %>%
        left_join(tinfo, by='target') %>% 
        mutate(pct=round(n_match/n_target*100, 2)) %>%
        arrange(desc(pct), desc(n_match)) %>%
        distinct(gene_id, .keep_all = T) %>%
        select(gene_id, target, n_match, n_target, pct, s_mean, e_mean)
}

