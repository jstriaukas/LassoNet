\name{lasso.net.grid}
\alias{lasso.net.grid}
\title{
Estimates coefficients and connection signs over the grid values (\eqn{\lambda}1, \eqn{\lambda}2).
}
\description{
Function fits lasso and network regressions over the grid values of penalty coefficients. \eqn{\beta} updates are computed in R and C++. In addtion, function stores connection signs, number of iterations until convergence, convergence binary output and alternating signs in case divergence in connection matrix. 
}
\usage{lasso.net.grid(x,y ,beta.0,lambda1,lambda2,M1,m.iter,n.iter,iscpp=TRUE,tol,alt.num)}
\arguments{
\item{x}{\eqn{n \times p} input data matrix}
\item{y}{response vector or size  \eqn{n \times 1}}
\item{beta.0}{initial value for \eqn{\beta}. default - zero vector of size \eqn{n \times 1}}
\item{lambda1}{lasso penalty coefficient}
\item{lambda2}{network penalty coefficient}
\item{M1}{penalty matrix}
\item{m.iter}{maximum number of iterations for sign matrix updating. default - 100}
\item{n.iter}{maximum number of iterations for \eqn{\beta} updating. default - 1e5}
\item{iscpp}{binary choice for using cpp function in coordinate updates. 1 - use C++ (default), 0 - use R}
\item{tol}{convergence in \eqn{\beta} tolerance level. default - 1e-6}
\item{alt.num}{In case connection sign matrix alternates, alt.num last iterataions are stored. default - 12 }
}
\details{Function loops through (using two for loops and one while loop) grid values of \eqn{\lambda}1 and \eqn{\lambda}2 until either sign matrix convergences or maximum (m.iter) number of iterations is reached. Convergence in sign matrix is checked by evaluating whether two consecutive matrices are equal. Warm starts are stored for each iterator. The first warm start is stored once coordinate updates converge given initial value for sign matrix (i.e. first iteration in sign matrix). Warm start for a new value of \eqn{\lambda}2 is stored when coordinate updates converge given initial value for sign matrix (i.e. first iteration in sign matrix, but last iteration in \eqn{\lambda}1). In addition, function checks whether signs alternate. Once the algorithm is close to the maximum number of iterations (m.iter) in sign matrix, function stores two previous sign matrices and checks whether they
alternate (wheter sign matrix equals to the previous two steps estimated sign matrix and does not equal to one step before estimated sign matrix). Further convergence checks are done for \eqn{\beta} values when looping though sign matrix. If at least one \eqn{\beta} updating procedure did not converge, \code{false} is stored in a convergece.in.grid output matrix while otherwise \code{true} is stored.}
\value{
\item{beta}{matrix of \eqn{\beta} coefficients. Columns denote different \eqn{\lambda}1
coefficients, rows - \eqn{\lambda}2 coefficients}
\item{mse}{Mean squared error value}
\item{M}{Array of connection signs. \eqn{M[,,i,j]} is the connection sign matrix for j-th \eqn{\lambda}1 value and i-th \eqn{\lambda}2 value}
\item{iterations}{matrix with stored number of steps for sign matrix to converge}
\item{update.steps}{matrix with stored number of steps for \eqn{\beta} updates to converge. (only stores the last values from connection signs iterations)}
\item{convergence.in.M}{matrix with stored values for convergence in sign matrix}
\item{convergence.in.grid}{matrix with stored values for convergence in \eqn{\beta} coefficients. If at least one \eqn{\beta} did not converge in sign matrix iterations, 0 (false) is stored, otherwise 1 (true)}
\item{xi.conv}{array with stored values about changes in connection signs (each iteration)}
\item{beta.alt}{vector for \eqn{\beta} value in the case when M1 matrix alternates, zero otherwise}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{
Weber, M., Schumacher, M., & Binder, H. (2014) \emph{Regularized Regression Incorporating Network Information: Simultaneous Estimation of Covariate Coefficients and Connection Signs}  

Striaukas, J. & Weber, M. (2018) \emph{Network Constrained Covariate Coefficient and Connection Sign Estimation}
}

\examples{
p=200
n=100
beta.0=array(1,c(p,1))
x=matrix(rnorm(n*p),n,p)
y=rnorm(n,mean=0,sd=1)
lambda1=c(0,1)
lambda2=c(0,1)
M1=diag(p)
lasso.net.grid(x, y, beta.0, lambda1, lambda2, M1)
}
