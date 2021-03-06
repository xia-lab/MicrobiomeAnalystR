% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_io.R
\name{.read.microbiomeanalyst.lib.rds}
\alias{.read.microbiomeanalyst.lib.rds}
\title{Read RDS files from the internet}
\usage{
.read.microbiomeanalyst.lib.rds(filenm, sub.dir = NULL, ref = NA)
}
\arguments{
\item{filenm}{Input the name of the file to download.}

\item{ref}{Default set to "NA".}

\item{opt}{Default set to "none".}
}
\description{
Function downloads the required file and reads it only if not already in working directory.
Need to specify the file URL and the destfile.
}
