% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SignalMatrix.R
\docType{methods}
\name{getCorrRegions}
\alias{getCorrRegions}
\title{# Filter rows in SignalMatrix by the correlation with their nearby genes}
\usage{
getCorrRegions(
  SignalMatrix,
  signal.granges,
  gene.matrix,
  corr.method = "spearman",
  species = "hg38",
  rule = "basalPlusExt",
  plot = FALSE,
  cores = 1
)
}
\arguments{
\item{signal.granges}{a path to a gtf file}

\item{gene.matrix}{a numeric matrix of gene expression with the same dimensions as `signal`}

\item{corr.method}{a character string indicating which correlation coefficient
is to be used for the test. One of "pearson", "kendall", or
"spearman", can be abbreviated.}

\item{species}{species name supported by rGREAT such as "hg38", "hg19", "mm10", "mm9" etc.}

\item{rule}{a rule how to associate genomic regions to genes supported by RGREAT, such as
"basalPlusExt", "twoClosest", "oneClosest"}

\item{plot}{a boolean value indicating whether rGREAT::plotRegionGeneAssociationGraphs()
will generates plots}

\item{cores}{a number of cores}

\item{sig.mtx}{a SignalMatrix object or a numeric matrix}
}
\value{
returns a \code{Model} object
}
\description{
# Filter rows in SignalMatrix by the correlation with their nearby genes
}
\details{

}
\examples{

}
