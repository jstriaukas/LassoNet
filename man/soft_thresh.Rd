\name{soft.thresh}
\alias{soft.thresh}
\title{
Soft thresholding operator.
}
\description{
Evaluates each element of \eqn{\beta} using soft thresholding operator.
}
\usage{soft.thresh(x, kappa)}
\arguments{
\item{x}{\eqn{\beta} coordinate}
\item{kappa}{\eqn{\kappa} value in general or \eqn{\lambda_1} for covariance updating.}
}

\details{Uses standard soft thresholding definition: \eqn{S(x,\kappa) = sign(x) (|x| - \kappa)_{+}}}
\value{
\item{x}{value after applying soft thresholding operator}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}

\examples{
kappa<-0.2
x<-0.7
soft.thresh(x, kappa)
}
