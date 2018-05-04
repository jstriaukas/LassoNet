\name{LassoNet-package}
\alias{LassoNet-package}
\docType{package}
\title{
LassoNet: package for 3CoSE algorithm. 

}
\description{
LassoNet R package estimates a penalized regression model slope coefficients together with connection signs between covariates using the 3CoSE algorithm described in the paper Weber et. al(2014) and briefly discussed in Striaukas and Weber (2018). Regression coefficients are penalized both in \eqn{\ell_1} and \eqn{\ell_2} norms, where the latter uses a network penalty matrix. The algorithm uses warm starts for the \eqn{\beta} vector and connection sign matrix \eqn{M}. The main function of the package is lasso.net.grid, see example below. Also, the package has a separate file with functions to fully replicate the results.}
\details{
\tabular{ll}{
Package: \tab LassoNet\cr
Type: \tab Package\cr
Version: \tab 0.8.3\cr
Date: \tab 2018-04-27 \cr
License: \tab Open source \cr
}
}

\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{
Weber, M., Schumacher, M., & Binder, H. (2014) \emph{Regularized Regression Incorporating Network Information: Simultaneous Estimation of Covariate Coefficients and Connection Signs}  

Striaukas, J. & Weber, M. (2018) \emph{Network Constrained Covariate Coefficient and Connection Sign Estimation}
}
\seealso{
\code{\link[Rcpp:Rcpp-package]{Rcpp}}, \code{\link[glmnet:glmnet-package]{glmnet}}
}
