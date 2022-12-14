% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualise.R
\docType{methods}
\name{plotMotifsActivity}
\alias{plotMotifsActivity}
\title{# Plot multiple motifs activities}
\usage{
plotMotifsActivity(
  motif.names,
  coef,
  title = motif.names,
  colour.line = "black",
  colour.xaxis = "black",
  ylab = "Motif activity",
  legend.position = "top",
  ...
)
}
\arguments{
\item{coef}{a matrix of coefficients. In columns are samples, in rows motifs.}

\item{title}{a character vector indicating a title of the plot}

\item{colour.line}{a character indicating a colour of a line}

\item{colour.xaxis}{a character indicating a colour of labels on X axis}

\item{legend.position}{a character indicating a position of legend, such as 
"left", "top", "right", or "bottom".}

\item{motif.name}{a character indicating a name of a DNA motif}
}
\value{
returns a \code{gtable} object from a ggplot R package
}
\description{
# Plot multiple motifs activities
}
