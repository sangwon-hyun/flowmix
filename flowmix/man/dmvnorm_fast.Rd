% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eigen.R
\name{dmvnorm_fast}
\alias{dmvnorm_fast}
\title{From eigendecomposition of the sigmas, calculate the same thing as
\code{mvtnorm::dmvnorm()}.}
\usage{
dmvnorm_fast(y, mu, sigma_eig)
}
\arguments{
\item{y}{Multivariate data.}

\item{mu}{Mean vector.}

\item{sigma_eig}{Result of eigendecomposition of sigma (using
\code{eigendecomp_sigma(sigma)} for covariance matrix \code{sigma}).}
}
\value{
Density vector.
}
\description{
From eigendecomposition of the sigmas, calculate the same thing as
\code{mvtnorm::dmvnorm()}.
}
