% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.RSEM.R
\name{RSEM_read_logFCs}
\alias{RSEM_read_logFCs}
\title{RSEM: read logFC data}
\usage{
RSEM_read_logFCs(files, labs = NULL)
}
\arguments{
\item{files}{file paths to "compare" files generated by \code{\link{RSEM_pipeline}}}

\item{labs}{character vector.}
}
\value{
tibble
}
\description{
Read and combine logFC.
}
\details{
Combie logFC for downstream analysis.
}
\author{
ZG Zhao
}
