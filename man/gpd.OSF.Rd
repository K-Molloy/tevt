% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gposf.R
\name{gpd.OSF}
\alias{gpd.OSF}
\title{Handler function for a series of Bias-based procedures for Choosing the Optimal Sample Fraction}
\usage{
gpd.OSF(data, method = "DK", ...)
}
\arguments{
\item{data}{vector of sample data}

\item{method}{selected method \code{default} is set to \code{DK}}

\item{...}{extra parameters for respective functions}
}
\value{
\item{second.order.par }{gives an estimation of the second order parameter \code{rho}.}
\item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
\item{threshold}{the corresponding threshold}
\item{tail.index}{the corresponding tail}
}
\description{
Handler function for a series of Bias-based procedures for Choosing the Optimal Sample Fraction
}
\details{
Modifying the method parameter allows access to 4 similar bias based procedures
}
\examples{
data(danish)
gpd.OSF(danish)
}
