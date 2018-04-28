#' Computes decomposition elements, p.g. 14 of Weber et. al (2014).
#' 
#' @param x Input data matrix of size \eqn{n \times p}, n - number of observations, p - number of covariates
#' @param y Response vector or size \eqn{n \times 1}
#' @param M Penalty matrix
#' @return Bx array of \eqn{B^{ij}_X} stored matrices. \eqn{Bx[,,k]} are the k-th combination of i and j non zero entry in the penalty matrix M
#' @return By array of \eqn{B^{ij}_y} stored matrices. \eqn{By[,k]} are the k-th combination of i and j non zero entry in the penalty matrix M
#' @examples
#' p<-200
#' n<-100
#' x<-matrix(rnorm(n*p),n,p)
#' y<-rnorm(n,mean=0,sd=1)
#' M<-diag(p)
#' get.BxBy(x, y, M)
get.BxBy   <- function(x , y ,M){
  # check only lower triangular values, since M is symmetric
  M[upper.tri(M)] <- 0
  non.zero <- which(M!=0,arr.ind = T)
  # in case diagonal has a non zero entry, we exlude that
  non.diag <- non.zero[,1]!=non.zero[,2]
  # store indicies of lower triagular signs being not zero
  e.con    <- non.zero[non.diag,]
  # number of iterations needed 
  k        <- nrow(e.con)
  # number of covariates
  p        <- ncol(x)
  # initializing
  By       <- array(0, c(2,k))
  Bx       <- array(0, c(2,p-2,k))
  if (k!=0){ # if there are off diagonal elements (i.e. if we are not at the ridge case)
    # loop through all indicies
    for (d in 1:k){
      elem    <- e.con[d,]
      
      x.tilde <- cbind(x[,elem[1]],x[,elem[2]])
      x.exc   <- cbind(x[,-c(elem[1],elem[2])])
      
      Byij    <- fastols(y    , x.tilde)
      Bxij    <- fastols(x.exc, x.tilde)
      
      By[,d]  <- Byij
      Bx[,,d] <- Bxij
    }
  } else{  
    Bx <- 0
    By <- 0
  }
  return(list(Bx = Bx, By = By))
}


get.signs.M <- function(MAT) {
  vec.out <- cbind(which(MAT!=0, arr.ind = T),as.matrix(MAT[which(MAT!=0, arr.ind = T)]))
  return(vec.out)
}


# update xi (actually used)
get.xi <- function(Bx ,By ,beta ,xi ,M){
  
  # check only lower triangular values, since M is symmetric
  M[upper.tri(M)] <- 0
  non.zero <- which(M!=0,arr.ind = T)
  # in case diagonal has a non zero entry, we exlude that
  non.diag <- non.zero[,1]!=non.zero[,2]
  # store indicies of lower triagular signs being not zero
  e.con    <- non.zero[non.diag,]
  # number of iterations needed 
  dd       <- nrow(e.con)
  if (dd!=0){
    for (i in 1:dd){
      beta.exc   <- beta[-c(e.con[i,1],e.con[i,2])]
      beta.tilde <- By[,i] - Bx[,,i]%*%beta.exc
      # updating only lower triangular part
      xi[e.con[i,1],e.con[i,2]] <- sign(beta.tilde[1])*sign(beta.tilde[2])
      
    }
    # updating upper triagular part
    xi[upper.tri(xi)] <- t(xi)[upper.tri(xi)]
    diag(xi) <- 1
  } else {
    xi <- M
  }
  return(xi)
}

# update M with a given initial M1 and estimated xi
matrix.M.update <- function(M,xi){
  init       <- -(xi)*abs(M) # update M1 matrix
  diag(init) <- abs(diag(M)) # store diagonal values as in initial penalty matrix
  return(init)
}


# Function for soft-thresholding operator
soft.thresh <- function(x ,kappa){
  # initializing
  
  # changing according to rule
  if (abs(x)>kappa){
    x <- sign(x)*(abs(x)-kappa)
  } else{
    x <- 0
  }
  
  return(x)
}


###### coefficient updates ##########



beta.update.net <- function(x ,y ,beta=array(0,length(y)) ,lambda1 ,lambda2 ,M1 ,n.iter = 1e5 ,iscpp=TRUE ,tol=1e-6 ){
  
  beta <- as.matrix(beta)
  x    <- as.matrix(x)
  M1    <- as.matrix(M1)
  
  # initializing the variables
  iii    <- length(beta)
  beta.previous           <- numeric(length(beta))
  k <- 1
  conv <- FALSE
  xy <- crossprod(x,y)
  xx <- crossprod(x)
  
  if (iscpp==TRUE){ # using cpp
    # access cpp function
    est   <- betanew_lasso_cpp(xx, xy, beta, M1, y, lambda1, lambda2, n.iter, tol)
    # store beta
    beta  <- est$beta
    # store steps until beta convergence
    k    <- est$steps
    # store value for convergence
    if(k < n.iter){
      conv <- TRUE
    }
    k <- k+1 # cpp starts at 0, thus the output gives one step less than R
  } else { # using R 
    
    while ((conv == FALSE) && (k-1 < n.iter) ) {
      # store previous beta
      beta.previous   <-  beta
      # update each at a time each coordinates
      for (j in 1:iii){
        #  sum_i=1^n [ 2*{n*beta.tilde(j) + <xj,y> - sum_k=1^p [ <xj,xk> beta.tilde(k) ]} - 2*lambda2 sum_j!v [ M(jv) beta.tilde(v) ]
        numerator   <- 2*((length(y)*beta[j] + xy[j] - xx[,j]%*%beta)) - 2*lambda2*M1[j,-j]%*%beta[-j]
        #  2*n + 2*lambda2 M(jj)
        denominator <- 2*length(y)+2*lambda2*M1[j,j]
        # update beta(j) = numerator/denominator
        beta[j] <- soft.thresh(numerator,lambda1)/(denominator)
      }
      #update beta
      conv <- sum(abs(beta-beta.previous))/iii<tol
      k <- k+1
    }
  }
  return(list(beta = beta, convergence = conv, steps = k-1))
}



#### grid search function ####

lasso.net.grid <- function(x ,y ,beta.0 = array(0,length(y)) ,lambda1=c(0,1) ,lambda2=c(0,1) ,M1 ,m.iter = 100, n.iter = 1e5 , iscpp=TRUE, tol=1e-6,  alt.num=12){
  # initializing beta
  beta <- beta.0
  
  # define penalty grid space
  n1 <- length(lambda1)
  n2 <- length(lambda2)
  # convert to relevant data structure
  beta    <- as.matrix(beta)
  x       <- as.matrix(x)
  M1      <- as.matrix(M1)
  
  # initialize values
  signs   <- FALSE
  xi      <- sign(t(x)%*%x)
  M       <- matrix.M.update(M1, xi)
  k       <- 0
  alt.id  <- 0
  
  
  # storage
  beta.store  <- array(0,c(length(beta), n2, n1))   
  M.storej     <- NULL
  M.store      <- NULL
  k.store     <- array(0,c(n2,n1))
  n.store     <- array(0,c(n2,n1))
  conv.store  <- array(0,c(n2,n1)) 
  conv.store2 <- array(0,c(n2,n1))
  alt.store   <- array(0,c(dim(M1),alt.num))
  alt.beta    <- array(0,c(length(beta),alt.num))
  alt.xi.store<- array(0,c(dim(xi),alt.num))
  mse         <- array(0,c(n2,n1))
  
  alt.save    <- NULL
  xi.conv     <- NULL
  xi.conv1    <- NULL
  xi.conv2    <- NULL
  
  
  # storage for warm starts
  beta.l1w <- numeric(length(beta)) # warm start for lambda1 
  beta.l2w <- numeric(length(beta)) # warm start for lambda2
  
  # store matrices Bx and By for xi update before loops 
  BxBy <- get.BxBy(x,y,M1)
  Bx   <- BxBy$Bx
  By   <- BxBy$By
  
  
  
  for (j in 1:n1) {
    # loop through lambda1
    
    cat("Looping through lambda1 - ", j, "\n") 
    for (i in 1:n2) {
      
      # loop through lambda2
      cat("Looping through lambda2 - ", i, "\n") 
      # store check value for beta convergence given updated M 
      betaconv <- FALSE
      # while loop for xi convegence
      while ((signs != TRUE)  && (k < m.iter) && (betaconv != TRUE)) {
        
        # estimates lasso with network given lambda1, lambda2 and M matrix
        est <- beta.update.net(x,y,beta,lambda1[j],lambda2[i],M,n.iter,iscpp)
        
        # update xi with "covariance" type updates. sec 2.2.
        xi.u <- get.xi(Bx,By,est$beta,xi,M1)
        
        # store coordinate descent convergence output
        if(est$convergence == TRUE) {
          betaconv <- TRUE
        } 
        
        # check for alternating signs if iterator is close to maximum number of iteration
        if(k >= m.iter-alt.num){
          id <- m.iter-k-alt.num+1
          alt.xi.store[,,id] <- xi
          alt.beta[,id]     <- est$beta
        } 
        if(k == m.iter) 
        {
          alt.id                   <- alt.id + 1 
          alt.save$beta[[alt.id]]  <- alt.beta
          alt.save$xi[[alt.id]]    <- alt.xi.store
          alt.save$lams[[alt.id]]  <- c(lambda1[j],lambda2[i])
        }
        
        # update iterations number
        k <- k + 1
        # check xi matrix with all previous matrices
        xi.conv[[k]] <- cbind(which(xi.u!=xi, arr.ind = T),as.matrix(xi.u[which(xi.u!=xi, arr.ind = T)]))
        # check for convergence or stop once interations reaches n.iter
        signs <- isTRUE(all.equal(xi.u,xi))
        
        # update xi
        xi <- xi.u
        # storing warm starts for beta lambda1 (the first converged beta)
        if (k == 1 && j == 1){
          beta.l1w <- est$beta
        }
        # updating beta for lambda2 warm start
        beta <- est$beta
        if(betaconv == FALSE) { 
          cat("coordinate descent did not converge", beta, "/n")
        }
      } # Iterations through sign Matrix finishes here. [either converged or signs alternate]
      
      
      xi.conv1[[i]] <- xi.conv
      
      k.store[i,j]     <- k 
      n.store[i,j]     <- est$steps
      conv.store[i,j]  <- signs
      
      if(betaconv == FALSE) { 
        conv.store2[i,j] <- FALSE 
      }
      else {
        conv.store2[i,j] <- betaconv 
      }
      
      k     <- 0
      signs <- FALSE
      
      
      
      beta.store[,i,j]  <- beta
      # compute MSE 
      mse[i,j]  <- mean((y - x%*%beta)^2)
      # store connections
      M.storej[[i]]   <- get.signs.M(M)
      # upate the penalty matrix 
      M   <- matrix.M.update(M1,xi)
      
    } # end of for loop n2  (through lambda2)
    beta <- beta.l1w # warm start for lasso
    M.store[[j]] <- M.storej
    # store xi check
    xi.conv2[[j]] <- xi.conv1
    # update beta for lambda1 warm start
    beta <- beta.l2w
  }   # end of for loop n1  (through lambda1)
  
  return(list(beta = beta.store, mse = mse, signMat = M.store, iterations = k.store, update.steps = n.store, convergence.in.M = conv.store, convergence.in.grid = conv.store2, xi.conv = xi.conv2, alt.info = alt.save))
}


lasso.net.fixed <- function(x ,y ,beta.0 = array(0,length(y)) ,lambda1=c(0,1) ,lambda2=c(0,1) ,M1 ,n.iter = 1e5 ,iscpp=TRUE ,tol=1e-6) {
  # initializing beta
  beta <- beta.0
  
  n1 <- length(lambda1)
  n2 <- length(lambda2)
  # convert to relevant data structure
  beta   <- as.matrix(beta)
  x      <- as.matrix(x)
  M1     <- as.matrix(M1)
  
  # storage
  beta.store  <- array(0,c(length(beta), n2, n1))   
  mse         <- array(0,c(n2,n1))
  n.store     <- array(0,c(n2,n1))
  conv.store2 <- array(0,c(n2,n1))
  # storage for warm starts
  beta.l1w <- numeric(length(beta)) # warm start for lambda1 
  beta.l2w <- numeric(length(beta)) # warm start for lambda2
  g        <- 1
  
  for (j in 1:n1) {
    # loop through lambda1
    
    cat("Looping through lambda1 - ", j, "\n") 
    for (i in 1:n2) {
      
      # loop through lambda2
      cat("Looping through lambda2 - ", i, "\n") 
      
      betaconv <- FALSE
      
      
      while(betaconv!=TRUE){
        est <- beta.update.net(x,y,beta,lambda1[j],lambda2[i],M1,n.iter,iscpp)
        # store coordinate descent convergence output
        if(est$convergence == TRUE) {
          betaconv <- TRUE
        } 
        beta <- est$beta
        # store beta
        beta.store[,i,j]  <- beta
        # compute MSE 
        mse[i,j]  <- mean((y - x%*%beta)^2)
        # store number of steps until convergence
        n.store[i,j]     <- est$steps
        if (j == 1){
          beta.l1w <- est$beta
        }
        if(betaconv == FALSE) { 
          cat("coordinate descent did not converge", beta, "/n")
        }
        if(betaconv == FALSE) { 
          conv.store2[i,j] <- FALSE 
        }
        else {
          conv.store2[i,j] <- betaconv 
        }  
      } 
    } # end of loop (n2)
    beta  <- beta.l1w # warm start for lasso
  }  # end of loop (n1)
  
  return(list(beta = beta.store, mse = mse, update.steps = n.store, convergence.in.grid = conv.store2))
  
}


# draw data

drawdata <- function(n,beta){
  p <- length(beta)
  x <- matrix(rnorm(n*p),n,p)
  y <- x%*%beta + rnorm(n, mean = 0, sd = sqrt(sum(beta^2)/4))
  return(list(y=y,x=x))
}


# fast ols (file fast_ols shows with simulated data that this is the fastest way in R)
fastols <- function (y ,x) { 
  XtX <- crossprod(x) 
  Xty <- crossprod(x, y) 
  return(solve(XtX, Xty)) 
}


mat.to.laplacian <- function(M1, type = "normalized"){
  
  degree.vec <- apply(M1,2,sum)
  
  if(type=="normalized"){
    du.dv.mat <- sqrt(as.matrix(degree.vec) %*% t(as.matrix(degree.vec)))
    du.dv.mat[which(du.dv.mat == 0)] <- 1
    L <- -1 * M1 / du.dv.mat
    L <- L + diag(sign(degree.vec))
  } 
  if(type=="combinatorial"){
    L <- - M1 + diag(degree.vec)
  }
  return(L)
}



