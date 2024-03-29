% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmp_utils.R
\name{PerformDEAnalyse}
\alias{PerformDEAnalyse}
\title{Perform differential analysis}
\usage{
PerformDEAnalyse(
  mbSetObj,
  taxalvl = "Genus",
  netType = "gem",
  overlay,
  initDE = 1,
  analysisVar = "CLASS",
  adjustedVar,
  alg = "limma",
  plvl = 0.05,
  fc.lvl = 1,
  selected = "NA",
  nonpar = FALSE
)
}
\arguments{
\item{mbSetObj}{Input the name of the mbSetObj.}

\item{taxalvl}{Character, input taxonomy level}

\item{metadata}{Character, input the name of the experimental factor
to group the samples.}
}
\description{
This functions performs alpha diversity.
}
\author{
Jeff Xia \email{jeff.xia@mcgill.ca}
McGill University, Canada
License: GNU GPL (>= 2)
}
