\name{matrix.M.update}
\alias{matrix.M.update}
\title{
Updates connection sign matrix.
}
\description{
Function updates M matrix using relation \eqn{(M)_{ij} = -\hat{\xi}_{ij} |(M_1)|_{ij}}.
}
\usage{matrix.M.update(M, xi)}
\arguments{
\item{M}{penalty matrix}
\item{xi}{estimated \eqn{\hat{\xi}_{ij}} matrix}
}
\details{Updating equation is at Weber et. al (2014), p.g. 9, 2nd step of algorithm}
\value{
\item{M}{updated correct M matrix}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{
Weber, M., Schumacher, M., & Binder, H. (2014) \emph{Regularized Regression Incorporating Network Information: Simultaneous Estimation of Covariate Coefficients and Connection Signs}  

Striaukas, J. & Weber, M. (2018) \emph{Network Constrained Covariate Coefficient and Connection Sign Estimation}
}
\examples{
p<-100
M<-diag(p)
xi<-matrix(rnorm(p*p), p, p)
matrix.M.update(M,xi)
}
