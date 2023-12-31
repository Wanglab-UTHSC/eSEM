% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pspMotifAdjacency.R
\name{pspMotifAdjacency}
\alias{pspMotifAdjacency}
\title{pspMotifAdjacency}
\usage{
pspMotifAdjacency(sp, ep, databasePSP, motif.ref)
}
\arguments{
\item{sp}{Substrate species; ("mouse", "rat", "human").}

\item{ep}{Enzyme species; ("mouse", "rat", "human").}

\item{databasePSP}{Default database or customized database input.}

\item{motif.ref}{Added enzyme-subatrate relationships from motif discovery. Same as KSEM -motif parameter.}
}
\description{
This script is used to build enzyme-substrate adjacency matrix with PSP database + Motif.
}
