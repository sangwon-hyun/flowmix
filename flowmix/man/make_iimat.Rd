% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.R
\name{make_iimat}
\alias{make_iimat}
\title{Indices for the cross validation jobs.}
\usage{
make_iimat(cv_gridsize, nfold, nrep)
}
\arguments{
\item{cv_gridsize}{CV grid size.}

\item{nfold}{Number of of CV folds.}

\item{nrep}{Number of repetitions.}
}
\value{
Integer matrix.
}
\description{
The resulting iimat looks like this:
}
\details{
ind ialpha ibeta ifold irep
 55      6     1     2    1
 56      7     1     2    1
 57      1     2     2    1
 58      2     2     2    1
 59      3     2     2    1
 60      4     2     2    1
}
