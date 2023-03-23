#' @name GOFULL
#' @aliases KEGGFULL
#' @title General GO and KEGG annotation data
#' @description `GOFULL` is parsed from 'http://purl.obolibrary.org/obo/go.obo'. `KEGGFULL` is from 'https://www.genome.jp/kegg-bin/get_htext?ko00001'.
#' @details
#' KEGGFULL and GOFULL are used for enrichment analysis in function \code{\link{statKEGG}} and \code{\link{statGO}}.
#' All data are tibble objects. Run examples for more details of the object.
#' @author ZG Zhao
#' @examples
#' data(GOFULL, package='vizRNAseq')
#' str(GOFULL)
#' data(KEGGFULL, package='vizRNAseq')
#' str(KEGGFULL)
NULL

#' @name ATHGENEANNO
#' @aliases ATHGENGO
#' @title Arabidopsis: gene annotation data
#' @description Both ATHGENEANNO and ATHGENGO are parsed from 'https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/gene_association.tair.gz'.
#' @details
#' Both `ATHGENEANNO` and `ATHGENGO` have `gene_id` column.
#' The `desc` column of `ATHGENEANNO` data contains a bunch of informations seperated by `|`. You can extract what you want.
#' All data are tibble objects. Run examples for more details of the object.
#' @author ZG Zhao
#' @examples
#' library(vizRNAseq)
#' str(ATHGENEANNO)
#' str(ATHGENEGO)
NULL
