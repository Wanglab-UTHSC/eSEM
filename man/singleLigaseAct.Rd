% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singleLigaseAct.R
\name{singleLigaseAct}
\alias{singleLigaseAct}
\title{singleLigaseAct}
\usage{
singleLigaseAct(ligase, input, adj, cor.off, kmo.off)
}
\arguments{
\item{ligase}{Single ligase name.}

\item{input}{Normalized and transformed ubiquiti-proteomics data.}

\item{adj}{Adjacency matrix.}

\item{cor.off}{Set up correlation cutoff value 0-1 to remove high collinear variables. Default is 0.95.}

\item{kmo.off}{Set up KMO cutoff value 0-1. Default is 0.}
}
\description{
Get single ligase activity.
}
