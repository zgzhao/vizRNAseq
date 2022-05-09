#' theme_sci
#'
#' ggplot2 theme_sci, a clear theme for figures of scientific papers.
#' @title theme_sci ggplot2 theme.
#' @param ... args passed to theme_bw().
#' @return ggplot2 theme set
#' @author ZG Zhao
#' @export
theme_sci <- function(base_size = 12, base_family = "Helvetica", ...) {
    theme_bw(base_size, base_family, ...) +
    theme(panel.grid=element_blank(),    #去掉网格线
          legend.title=element_blank(),  #去掉图例title
          legend.key=element_rect(color="white",fill="white")
      )
}

#' @export
style_xtext <- function(angle=30, hjust=1, vjust=1, ...) {
    theme(axis.text.x=element_text(angle=angle, hjust=hjust, vjust=vjust, ...))
}

#' @export
style_ytext <- function(...) {
    theme(axis.text.y=element_text(...))
}

#' @export
style_default_bar <- function(horiz=FALSE, cex=0.1) {
    if(! horiz) {
        scale_y_continuous(expand = expansion(mult=c(0, cex)))
    } else {
        scale_x_continuous(expand = expansion(mult=c(0, cex)))
    }
}
