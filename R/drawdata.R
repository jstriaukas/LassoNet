#' Draws random data
#' 
#' @param n number of observations
#' @param beta slope parameter vector
#' @return y response vector
#' @return x design matrix
drawdata <- function(n,beta){
  p <- length(beta)
  x <- matrix(rnorm(n*p),n,p)
  y <- x%*%beta + rnorm(n, mean = 0, sd = sqrt(sum(beta^2)/4))
  return(list(y=y,x=x))
}