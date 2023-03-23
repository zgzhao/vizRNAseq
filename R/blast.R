mseq_len <- function(ss, ee) {
    xx <- IRanges(start=ss, end = ee)
    width(IRanges::reduce(xx))
}

.reduceRead <- function(fx, ids) {
    ftmp1 <- tempfile()
    ftmp2 <- tempfile()
    on.exit({
        if(file.exists(ftmp1)) file.remove(ftmp1)
        if(file.exists(ftmp2)) file.remove(ftmp2)
    })
    writeLines(as.character(unique(ids)), ftmp1)
    pperl <- file.path(path.package("vizRNAseq"), "Perl", "reduce-table.pl")
    system(paste("perl", pperl, ftmp1, fx, ftmp2))
    read.table(ftmp2)
}

#' Parsing blast results for best matches. Perl programs included: make sure you have Perl installed in your system.
#'
#' Best results are those: (1) have longest matches in query sequence; (2) larger bit_score if matches are equal.
#' @title Get best matches from blast results.
#' @param blast_result character, file path to BLAST output file in fmt6 format.
#' @return tibble
#' @author ZG Zhao
#' @export
bestBlast <- function(blast_result) {
    xfile <- tempfile()
    on.exit({
        if(file.exists(xfile)) file.remove(xfile)
    })
    pperl <- file.path(path.package("vizRNAseq"), "Perl", "best-blast.pl")
    system(paste("perl", pperl, blast_result, xfile))
    dtx <- read.table(xfile)
    colnames(dtx) <- c('gene_id', 'acc', 'len', 'score')
    return(dtx %>% as_tibble %>% arrange(gene_id))
}

#' Map sequences to GO or KO (KEGG) identifiers based on best blast results. Perl programs included: make sure you have Perl installed in your system.
#'
#' Besides a blast result file (fmt6 format), you should prepare a file mapping accession to GO or KO. If the accessions in the two files are of different types, an additional mapping file is required for name translation. To prepare the mapping files, You can download the data from NCBI (for GO mapping) or KEGG (for KO mapping) website. All files input must be tab-table without header. 
#' @title Pipeline: blast2go/blast2kegg
#' @aliases blast2kegg
#' @param file_fmt6 character, file path to NCBI-blast output in fmt6 format
#' @param file_acc_trans (optional) character, file path. Data should have two columns: orignal_id and final_id. Orignal ids are the target accessions in blast results and the final ids are those for GO or KO mapping.
#' @param file_acc2xo character, file path. Data should have two columns: accessions and GO/KO ids.
#' @param uniprot TRUE/FALSE. TRUE for uniprot accessions and FALSE (default) for ncbi accessions.
#' @param KEGG TRUE/FALSE (default). Switch to parsing KEGG results.
#' @param max_e double. Default value is 1e-5.
#' @param min_score numeric. Default value is 1.
#' @return tibble
#' @author ZG Zhao
#' @export
blast2go <- function(file_fmt6, file_acc_trans, file_acc2xo,
                     uniprot=FALSE, KEGG=FALSE) {
    ## step1: parse best blast
    cat("Parsing best blast ...\n")    
    ans <- bestBlast(file_fmt6) 
    
    if(uniprot) ans <- mutate(ans, acc=sub("^[^\\|]+\\|([^\\|]+)\\|.+$", "\\1", acc))
    ans <- mutate(ans, acc=sub("\\.[0-9]+$", "", acc))
    
    ## step2 (optional): translate accession ids
    if(! missing(file_acc_trans)) {
        cat("Translating accession ids ...\n")
        df_acc_trans <- .reduceRead(file_acc_trans, unique(ans$acc))
        colnames(df_acc_trans) <- c('acc', 'gid')
        ans <- ans %>%
            left_join(df_acc_trans, by='acc') %>%
            select(gene_id, gid) %>%
            rename(acc=gid)
    }
    
    ## step3: reduce and read annotation table
    cat("Parsing annotation ...\n")
    df_acc2xo <- .reduceRead(file_acc2xo, ans$acc)
    colnames(df_acc2xo) <- c('acc', 'GO')
    df_acc2xo <- df_acc2xo %>%
        group_by(acc) %>%
        summarise(GO=paste(unique(GO), collapse=" "))

    ## step4: merge data
    cat("Finalizing ...\n")
    ans <- ans %>%
        left_join(df_acc2xo, by='acc') %>%
        as_tibble %>%
        filter(! is.na(GO)) 
    
    ## step5: arrange and return results
    if(KEGG) ans <- rename(ans, ko=GO)
    ans
}

#' @export
blast2kegg <- function(...) blast2go(...)
