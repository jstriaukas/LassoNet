\name{get.BxBy}
\alias{get.BxBy}
\title{
Computes decomposition elements, p.g. 14 of Weber et. al (2014).
}
\description{
Function computes matrices \eqn{B^{ij}_X} and \eqn{B^{ij}_y} to speed up estimation of connection signs. These matrices are precomputed, stored. We store only indices that have non zero entries in penalty matrix M. 
}
\usage{get.BxBy(x, y, M)}
\arguments{
\item{x}{Input data matrix of size \eqn{n \times p}, n - number of observations, p - number of covariates}
\item{y}{y Response vector or size \eqn{n \times 1}}
\item{M}{Penalty matrix}
}


\details{Function calculates matrices all for i and j indices that have non zero
values in a given penalty matrix, since \eqn{\xi} matrix is updated only for
such indices.}

\value{
\item{Bx}{array of \eqn{B^{ij}_X} stored matrices. \eqn{Bx[,,k]} are the k-th combination of i and j non zero entry in the penalty matrix M}
\item{By}{array of \eqn{B^{ij}_y} stored matrices. \eqn{By[,k]} are the k-th combination of i and j non zero entry in the penalty matrix M}
}

\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{
Weber, M., Schumacher, M., & Binder, H. (2014) \emph{Regularized Regression Incorporating Network Information: Simultaneous Estimation of Covariate Coefficients and Connection Signs}  

Striaukas, J. & Weber, M. (2018) \emph{Network Constrained Covariate Coefficient and Connection Sign Estimation}
}

\examples{
p<-200
n<-100
x<-matrix(rnorm(n*p),n,p)
y<-rnorm(n,mean=0,sd=1)
M<-diag(p)
get.BxBy(x, y, M)
}
