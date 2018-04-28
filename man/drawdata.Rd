\name{drawdata}
\alias{drawdata}
\title{
Draw random data.
}
\description{
Function generates random data. 
}
\usage{drawdata(n, beta)}
\arguments{
\item{n}{number of observations}
\item{beta}{\eqn{\beta} coefficient vector}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}

\examples{
n<-10
beta<-rep(0.2,n)
data<-drawdata(n, beta)
}
