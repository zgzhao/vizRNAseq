% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.cuff.R
\name{mergeRead}
\alias{mergeRead}
\title{mergeRead function}
\usage{
mergeRead(
  files,
  d.col = "log2FC_adj",
  treats = NULL,
  sels = NULL,
  keep = TRUE,
  ...
)
}
\arguments{
\item{files}{filenames}

\item{d.col}{data column to merge, one column only.}

\item{treats}{Character vector of same length of filenames for names of treatments.}

\item{sels}{NULL. Character vector of gene_id to select.}

\item{keep}{TRUE/FALSE. keep the merge.by column in result data.frame.}

\item{...}{args passed to read.csv}
}
\value{
data.frame
}
\description{
Function for Cuff: read in and merge data.
}
\details{
Read data from cuffdiff result file and merge.
}
\author{
ZG Zhao
}
