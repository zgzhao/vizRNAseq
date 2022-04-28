#' parse blast result
#'
#' parse blast result.
#' @title Get best matches from blast results.
#' @param blast_result BLAST out put file in fmt6 format
#' @param blast_db BLAST db file (fasta) when you do blast 
#' @param prot TRUE/FALSE. Protein blast db (not query sequences).
#' @return tibble
#' @author ZG Zhao
#' @export
bestMatches <- function(blast_result, blast_db, prot=FALSE) {
    require(Biostrings)
    ss <- if(prot) readAAStringSet(blast_db) else readDNAStringSet(blast_db)
    xnames <- sub(' .+$', '', names(ss))
    tinfo <- tibble(target=xnames,
                    tlen=Biostrings::nchar(ss))

    dx <- read.table(blast_result) %>% as_tibble
    colnames(dx) <- c('gene_id', 'target', 'p_indentical',
                      'n_align', 'n_mismatch', 'gap_open',
                      'qstart', 'qend',
                      'sstart', 'send',
                      'evalue', 'bit_score')

    dx %>%
        mutate(gene_id = sub('_i[0-9]+$', '', gene_id),
               irs = map2(sstart, send, seq)) %>%
        group_by(gene_id, target) %>%
        summarise(n = length(unique(unlist(irs))),
                  x = sum(n_mismatch),
                  match = n - x) %>% 
        left_join(tinfo, by='target') %>% 
        mutate(pct=round(match/tlen*100, 2)) %>%
        arrange(desc(pct), desc(match)) %>%
        distinct(gene_id, .keep_all = T)
}

