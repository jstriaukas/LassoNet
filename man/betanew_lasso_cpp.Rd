\name{betanew_lasso_cpp}
\alias{betanew_lasso_cpp}
\title{
C++ subroutine that updates \eqn{\beta} coefficients.
}
\description{
This function updates \eqn{\beta} for given penalty parameters.
}
\usage{betanew_lasso_cpp(xx, xy, beta, M, y, Lambda1, Lambda2, iter, tol)}
\arguments{
\item{xx}{Bx matrix}
\item{xy}{By vector}
\item{beta}{initial value for \eqn{\beta}; default - zero vector of size \eqn{p \times 1}}
\item{M}{penalty matrix}
\item{y}{response vector or size \eqn{n \times 1}}
\item{Lambda1}{lasso penalty parameter}
\item{Lambda2}{network penalty parameter}
\item{iter}{maximum number of iterations for \eqn{\beta} step}
\item{tol}{convergence tolerance level}
}
\details{See beta.update.net}
\value{
\item{beta}{updated \eqn{\beta} vector}
\item{steps}{number of steps until convergence}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{
Weber, M., Schumacher, M., & Binder, H. (2014) \emph{Regularized Regression Incorporating Network Information: Simultaneous Estimation of Covariate Coefficients and Connection Signs}  

Weber, M., Striaukas, J., Schumacher, M., & Binder, H. (2018) \emph{Network Constrained Covariate Coefficient and Connection Sign Estimation}
}

\examples{
p<-200
n<-100
beta.0=array(1,c(p,1))
x<-matrix(rnorm(n*p),n,p)
y<-rnorm(n,mean=0,sd=1)
lambda1<-1
lambda2<-1
M1<-diag(p)
updates<-beta.update.net(x, y, beta.0, lambda1, lambda2, M1)
}
