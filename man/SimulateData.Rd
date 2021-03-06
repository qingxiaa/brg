% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulateData.R
\name{SimulateData}
\alias{SimulateData}
\title{This function is to simulate genes expression data for evaluating the performance of BRIDGE}
\usage{
SimulateData(
  J,
  N,
  G,
  M,
  TimeEff0,
  mu01,
  sigma01,
  sigma02,
  delta2,
  rho1,
  rho2,
  phi2
)
}
\arguments{
\item{J}{The number of genes}

\item{N}{The number of subjects at each time point}

\item{G}{The number of differential genes}

\item{M}{The number of bridging samples}

\item{TimeEff0}{Time effects across two time points for differential genes}

\item{mu01}{mean expression at time point 1}

\item{sigma01}{variance of gene expression at time point 1}

\item{sigma02}{variance of gene expression at time point 2}

\item{delta2}{multiplicative batch effects of batch 2}

\item{rho1}{correlation coefficient between different batch}

\item{rho2}{correlation coefficient between different time point}

\item{phi2}{additive batch effects of batch 2}
}
\value{
Function outputs a list including the following:
\describe{
\item{\code{Data}}{A simulated genomic measure matrix (observation x genes) with batch, time and subjectID information at first 3 columns }
\item{\code{bridgingSampleInd}}{The index of subjects that are bridging samples}
\item{\code{Geneind}}{The index of genes that are truely differential}
\item{\code{Obs_te}}{The vector of observed time effect size of truly differential genes}}
}
\description{
This function is to simulate genes expression data for evaluating the performance of BRIDGE
}
