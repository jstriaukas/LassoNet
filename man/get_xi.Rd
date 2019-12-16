\name{get.xi}
\alias{get.xi}
\title{
Updates the estimates of the connection signs by running mini OLS models.
\description{
Updates connection signs \eqn{\hat{\xi}}.
}
\usage{get.xi(Bx,By,beta,xi,M)}
\arguments{
\item{Bx}{Bx element }
\item{By}{By element}
\item{beta}{\eqn{\hat{\beta}} estimated value}
\item{xi}{\eqn{\hat{\xi}} matrix estimated at the previous step}
\item{M}{penalty matrix}
}
\details{}
\value{
\item{xi}{\eqn{\hat{\xi}} matrix}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{
Weber, M., Schumacher, M., & Binder, H. (2014). \emph{Regularized Regression Incorporating Network Information: Simultaneous Estimation of Covariate Coefficients and Connection Signs}  

Weber, M., Striaukas, J., Schumacher, M., & Binder, H. (2018) \emph{Network Constrained Covariate Coefficient and Connection Sign Estimation}
}
