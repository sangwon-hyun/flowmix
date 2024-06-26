% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin-data.R
\name{make_ybin_1d}
\alias{make_ybin_1d}
\title{The same as \code{make_ybin_3d()} but in 1d.}
\usage{
make_ybin_1d(counts, midpoints, names = NULL)
}
\arguments{
\item{counts}{d x d x d array containing counts. If this is NULL, then dummy
counts of -100 are added.}

\item{midpoints}{midpoints of each bin.}

\item{names}{Names of the data columns.}
}
\description{
The same as \code{make_ybin_3d()} but in 1d.
}
