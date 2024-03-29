% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gate.R
\name{gate_plot}
\alias{gate_plot}
\title{Using gated particles (from \code{gate()}) , create one plot of the top
levels.}
\usage{
gate_plot(
  tt,
  dims,
  y,
  counts,
  memberships,
  datetime = "",
  top_clusters,
  top_cluster_names = NULL
)
}
\arguments{
\item{tt}{Time point (one of \code{1:length(y_gated)}).}

\item{dims}{Which dimensions to plot.}

\item{y}{Particle data.}

\item{counts}{Particle multiplicities.}

\item{memberships}{Particle memberships.}

\item{datetime}{An object (string or datetime) containing the time and date,
to go in the title.}

\item{top_clusters}{Memberships (cluster numbers) you'd like to highlight.}

\item{top_cluster_names}{(optional) Labels for the clusters in
\code{top_clusters}.}
}
\value{
One ggplot object
}
\description{
Using gated particles (from \code{gate()}) , create one plot of the top
levels.
}
