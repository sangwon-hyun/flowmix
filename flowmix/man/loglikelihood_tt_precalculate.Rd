% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objective.R
\name{loglikelihood_tt_precalculate}
\alias{loglikelihood_tt_precalculate}
\title{First helper function: Calculates one particle's log likelihood using
precalculated data densities.}
\usage{
loglikelihood_tt_precalculate(
  ylist,
  tt,
  denslist_by_clust,
  prob,
  countslist = NULL,
  numclust
)
}
\arguments{
\item{ylist}{Response data.}

\item{tt}{Time point of interest.}

\item{denslist_by_clust}{Pre-calculated densities.}

\item{prob}{Population proportions.}

\item{countslist}{(Optional) Counts corresponding to \code{ylist}.}

\item{numclust}{Number of clusters.}
}
\value{
Log likelihood.
}
\description{
First helper function: Calculates one particle's log likelihood using
precalculated data densities.
}
