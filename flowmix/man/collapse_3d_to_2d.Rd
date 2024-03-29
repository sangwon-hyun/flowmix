% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-helpers.R
\name{collapse_3d_to_2d}
\alias{collapse_3d_to_2d}
\title{Collapses from a 3d cytogram to two dimensions. This is mainly used by
\code{scatterplot_2d()}.}
\usage{
collapse_3d_to_2d(y, counts, dims = 1:2)
}
\arguments{
\item{y}{3d cytogram (a 3-column matrix).}

\item{counts}{The multiplicity for each of the particles in \code{y}.}

\item{dims}{Two of \code{c(1:3)}.}
}
\value{
3-column matrix; first two columns are the dimensions in
  \code{dims}.
}
\description{
Collapses from a 3d cytogram to two dimensions. This is mainly used by
\code{scatterplot_2d()}.
}
