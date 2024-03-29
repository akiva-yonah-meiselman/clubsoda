% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesis_test.r
\name{css.test.lambdas}
\alias{css.test.lambdas}
\title{Eigenvalues for CSS Test}
\usage{
css.test.lambdas(xg, xtx, c.0)
}
\arguments{
\item{xg}{List of cluster-specific submatrices of design matrix}

\item{xtx}{Inverse of design matrix squared}

\item{c.0}{Linear hypothesis vector}
}
\value{
A vector of eigenvalues
}
\description{
This function gets the eigenvalues for a test based on
Carter, Schnepel, and Steigerwald (2017).
}
