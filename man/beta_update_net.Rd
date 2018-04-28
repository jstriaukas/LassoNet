\name{beta.update.net}
\alias{beta.update.net}
\title{
Updates \eqn{\beta} coefficients.
}
\description{
This function updates \eqn{\beta} coefficients for a given penalty parameters.
}
\usage{beta.update.net(x,y,beta,lambda1,lambda2,M1,n.iter,iscpp,tol)}
\arguments{
\item{x}{input data matrix of size \eqn{n \times p}, n - number of observations, p - number of covariates}
\item{y}{response vector or size \eqn{n \times 1}}
\item{beta}{initial value for \eqn{\beta}. default - zero vector of size \eqn{n \times 1}}
\item{lambda1}{lasso penalty coefficient}
\item{lambda2}{network penalty coefficient}
\item{M1}{penalty matrix}
\item{n.iter}{maximum number of iterations for \eqn{\beta} updating. default - 1e5}
\item{iscpp}{binary choice for using cpp function in coordinate updates. 1 - use C++ (default), 0 - use R.}
\item{tol}{convergence tolerance level. default - 1e-6}
}
\details{The function updates the coefficient vector \eqn{\beta} given the data and penalty parameters.
Convergence criterion is defined as: \eqn{\sum_{i=1}^p |\beta_{i,j} - \beta_{i,j-1}| \leq} tol.}
\value{
\item{beta}{updated \eqn{\beta} vector}
\item{convergence}{binary variable for convergence}
\item{steps}{number of steps until convergence}
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
beta.0=array(1,c(p,1))
x<-matrix(rnorm(n*p),n,p)
y<-rnorm(n,mean=0,sd=1)
lambda1<-1
lambda2<-1
M1<-diag(p)
updates<-beta.update.net(x, y, beta.0, lambda1, lambda2, M1)
}
