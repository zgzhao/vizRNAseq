#' Generate customized ontology file for KEGG network analysis in Cytoscape using BiNGO app.
#'
#' KEGG network analysis in Cytoscape using BiNGO app needs two additional files: a general ontology file and an species-specific annotation file. This function creates the former.
#' @title Function BiNGO_KEGG_ontology
#' @param out.file 
#' @return out.file written
#' @author ZG Zhao
#' @export
BiNGO_KEGG_ontology <- function(out.file){
    hds <- c("(curator=KEGG)(type=Metabolism)")
    kegx <- keggs
    kegx$description <- gsub(" *\\[.*$", "", kegx$description)
    kegx$ko <- sub("^K", "1", kegx$ko, ignore.case = T)
    anno <- apply(kegx, 1, FUN = function(x){
        isa <- setdiff(x[4:6], NA)
        if(length(isa) < 1) paste0(x[2], " = ", x[3])
        else paste0(x[2], " = ", x[3], " [isa: ", rev(isa)[1], " ]")
    })
    names(anno) <- NULL
    anno <- c(hds, anno)
    writeLines(anno, out.file)
    cat(paste("Ontology file for network analysis in Cytoscape BiNGO app:", out.file, "generated.\n"))
    invisible(anno)
}

#' Generate customized ontology file for KEGG network analysis in Cytoscape using BiNGO app.
#'
#' KEGG network analysis in Cytoscape using BiNGO app needs two additional files: an ontology file and an species-specific annotation file. Function "BiNGO_KEGG_annotation" is designed to generate a customized ontology file.
#' @title Function BiNGO_KEGG_annotation
#' @param df.anno gene to KO annotation data.frame, whose first column is gene id and second column is KEGG ko.
#' @param out.file path to save the final annotation file.
#' @param species.name character. Name of your species.
#' @return annotation file written
#' @author ZG Zhao
#' @export
BiNGO_KEGG_annotation <- function(df.anno, out.file, species.name){
    df.anno[[2]] <- toupper(df.anno[[2]])
    ss <- df.anno[[2]] %in% keggs$ko
    df.anno <- df.anno[ss, ]
    hline <- paste0("(species=", species.name, ")(type=Metabolism)(curator=KEGG)")
    df.anno[[2]] <- sub("^K", "1", df.anno[[2]], ignore.case = T)
    anno <- c(hline, paste(df.anno[[1]], df.anno[[2]], sep=" = "))
    writeLines(anno, out.file)
    cat(paste("Organism file for network analysis in Cytoscape BiNGO app:", out.file, "generated.\n"))
    invisible(anno)
}
