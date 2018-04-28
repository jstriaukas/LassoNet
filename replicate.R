#' 
#' run_simulations.single function runs a single scenario simulation 
#' 
run_simulations.single <- function(p,beta.x, beta.0, M, lambda.1, lambda.2, t, n, nfolds = 10,  tol = 1e-6, n.iter = 1e10, m.iter, mc = FALSE){
  
  if(!mc)OUTPUT <-lapply(1:t, main.loop.single, 
                         p,beta.x,
                         beta.0, M, lambda.1, lambda.2,
                         n, nfolds,  tol, n.iter, m.iter) 
  
  if(mc)OUTPUT <- sfLapply(1:t, main.loop.single, 
                           p,beta.x,
                           beta.0, M, lambda.1, lambda.2,  
                           n, nfolds,  tol, n.iter, m.iter) 
  return(OUTPUT)
}

simul.data.single      <- function(p, n, beta.x, q = 10) {
  
  X1  <- array(0,c(n,p))
  Y1 <- array(0,c(n,1))
  
  prop <- p/(q+1)
  # storage
  
  b1mat <- matrix(t(beta.x),q+1,prop)
  sig    <- 0.51
  for (jj in 1:n) {
    X1.temp      <- matrix(0,prop,q+1)
    TF          <- rnorm(prop,mean = 0, sd = 1)
    for (ii in 1:prop) {
      x.temp      <- rnorm(q,mean=rep(sign(b1mat[1,ii])*sign(b1mat[2:11,ii])*0.7*TF[ii],q), sd = sqrt(sig))
      X1.temp[ii,] <- c(TF[ii], x.temp)
    }
    X1[jj,]      <- as.vector(t(X1.temp))
  }
  # computing variance of each scenario:
  epsilon <- rnorm(n, mean = 0, sd = 1)
  sig.1   <- sqrt(sum((beta.x)^2)/4)
  Y1 <- X1%*%beta.x + sig.1*epsilon
  
  return(list(Y1 = Y1, X1 = X1))
}

main.loop.single <- function(dt,p,beta.x, beta.0, M, lambda.1, lambda.2, n, nfolds = 10, tol = 1e-6, n.iter = 1e10, m.iter = 2e1){
  
  
  
  #---params---#
  n1 <- length(lambda.1)
  n2 <- length(lambda.2)
  
  
  # storage 
  msfeCV.norm <- array(0,c(n2,n1,nfolds))
  msfeCV.comb <- array(0,c(n2,n1,nfolds))
  msfeCV.fix <- array(0,c(n2,n1,nfolds))
  msfeCV.lasso <- array(0,c(n2,n1,nfolds))
  
  OUTPUT        <- NULL
  MSFE.norm   <- NULL
  MSFE.comb  <- NULL
  MSFE.fix   <- NULL
  MSFE.lasso   <- NULL
  
  lam.comb.norm   <- NULL
  lam.comb.comb   <- NULL
  lam.comb.fix   <- NULL
  lam.comb.lasso   <- NULL
  
  std.cv.norm         <- NULL
  std.cv.comb         <- NULL
  std.cv.fix         <- NULL
  std.cv.lasso         <- NULL
  
  PRED.ERROR.NORM      <- NULL
  PRED.ERROR.COMB      <- NULL
  PRED.ERROR.FIX      <- NULL
  PRED.ERROR.LASSO      <- NULL
  
  PRED.ERROR.null <- NULL  
  PRED.ERROR.true <- NULL
  
  #--data train (cross validate) [i] and test [i+1]---#
  
  # i-th data set: training and validating:
  data.sim <- simul.data.single(p, n, beta.x) 
  X1       <- scale(data.sim$X1 , center = TRUE, scale = TRUE)
  Y1       <- scale(data.sim$Y1, center = TRUE, scale = FALSE)
  
  # i-th data set: testing and computing statistics:
  data.sim.test <- simul.data.single(p, n, beta.x) 
  X1.test       <- scale(data.sim.test$X1 , center = TRUE, scale = TRUE)
  Y1.test       <- scale(data.sim.test$Y1 , center = TRUE, scale = FALSE)
  
  
  
  if (nfolds < 2) 
    stop("nfolds must be bigger than 2; at least 3 would be good; nfolds=10 recommended")
  crossval.ind <- CVfold(n,K = nfolds, is.random = FALSE)
  
  CVOUTPUT <- cv.par(nfolds, dt, crossval.ind, X1, Y1, M, beta.0, lambda.1, lambda.2, m.iter, n.iter, tol)
  
  
  for(d in 1:nfolds){
    # store estimates of connections, convergence and info if signs alternate
    beta.norm   <- CVOUTPUT[[d]]$beta.norm
    beta.comb   <- CVOUTPUT[[d]]$beta.comb
    beta.fix    <- CVOUTPUT[[d]]$beta.fix
    beta.lasso  <- CVOUTPUT[[d]]$beta.lasso
    
    X1.val <- CVOUTPUT[[d]]$X1.val
    Y1.val <- CVOUTPUT[[d]]$Y1.val 
    
    
    for (jj in 1:n2) { # ii - lasso, jj - network
      for (ii in 1:n1){
        msfeCV.norm[jj,ii,d] <- mean((Y1.val - X1.val%*%as.numeric(beta.norm[,jj,ii]))^2)
        msfeCV.comb[jj,ii,d] <- mean((Y1.val - X1.val%*%as.numeric(beta.comb[,jj,ii]))^2)
        msfeCV.fix[jj,ii,d]  <- mean((Y1.val - X1.val%*%as.numeric(beta.fix[,jj,ii]))^2)
        msfeCV.lasso[1,ii,d] <- mean((Y1.val - X1.val%*%as.numeric(beta.lasso[,1,ii]))^2)
      }
    }
  }
  
  # average of CV error curves:
  OUTPUT$MSFE.norm   <- apply(msfeCV.norm, c(1,2), "mean")
  OUTPUT$MSFE.comb   <- apply(msfeCV.comb, c(1,2), "mean")
  OUTPUT$MSFE.fix    <- apply(msfeCV.fix, c(1,2), "mean")
  OUTPUT$MSFE.lasso  <- apply(msfeCV.lasso, c(1,2), "mean")
  
  OUTPUT$std.cv.norm   <- apply(msfeCV.norm, c(1,2), "sd")/sqrt(nfolds)
  OUTPUT$std.cv.comb   <- apply(msfeCV.comb, c(1,2), "sd")/sqrt(nfolds)
  OUTPUT$std.cv.fix    <- apply(msfeCV.fix, c(1,2), "sd")/sqrt(nfolds)
  OUTPUT$std.cv.lasso  <- apply(msfeCV.lasso, c(1,2), "sd")/sqrt(nfolds)
  
  
  lam.comb.norm  <- which(OUTPUT$MSFE.norm == min(OUTPUT$MSFE.norm), arr.ind = T)
  lam.comb.comb  <- which(OUTPUT$MSFE.comb == min(OUTPUT$MSFE.comb), arr.ind = T)
  lam.comb.fix   <- which(OUTPUT$MSFE.fix == min(OUTPUT$MSFE.fix), arr.ind = T)
  lam.comb.lasso <- which(OUTPUT$MSFE.lasso == min(OUTPUT$MSFE.lasso), arr.ind = T)
  
  OUTPUT$optim.lam.norm <- c(lambda.1[lam.comb.norm[2]], lambda.2[lam.comb.norm[1]])
  OUTPUT$optim.lam.comb <- c(lambda.1[lam.comb.comb[2]], lambda.2[lam.comb.comb[1]])
  OUTPUT$optim.lam.fix <- c(lambda.1[lam.comb.fix[2]], lambda.2[lam.comb.fix[1]])
  OUTPUT$optim.lam.lasso <- c(lambda.1[lam.comb.lasso[2]])
  
  
  
  cat('*** Norm Lap ***')
  M.int <- mat.to.laplacian(M, type = "normalized")
  OUT.NORM <- lasso.net.grid(X1,Y1,beta.0,OUTPUT$optim.lam.norm[1],OUTPUT$optim.lam.norm[2], M.int, m.iter = m.iter, n.iter = n.iter, iscpp = T, tol = tol, alt.num = 4 )
  OUTPUT$OUT.NORM <- OUT.NORM
  cat('*** Fixed Norm ***')
  M.int <- mat.to.laplacian(M, type = "normalized")  
  OUT.FIX <- lasso.net.fixed(X1,Y1,beta.0,OUTPUT$optim.lam.fix[1],OUTPUT$optim.lam.fix[2], M.int, n.iter = n.iter, iscpp = T, tol = tol )
  OUTPUT$OUT.FIX <- OUT.FIX
  cat('*** Comb Lap ***')
  M.int <- mat.to.laplacian(M, type = "combinatorial")
  OUT.COMB <- lasso.net.grid(X1,Y1,beta.0,OUTPUT$optim.lam.comb[1],OUTPUT$optim.lam.comb[2], M.int, m.iter = m.iter, n.iter = n.iter, iscpp = T, tol = tol, alt.num = 4 )
  OUTPUT$OUT.COMB <- OUT.COMB
  cat('*** Lasso ***')
  OUT.LASSO <- lasso.net.grid(X1,Y1,beta.0,OUTPUT$optim.lam.lasso[1],0, M, m.iter = m.iter, n.iter = n.iter, iscpp = T, tol = tol, alt.num = 4 )
  OUTPUT$OUT.LASSO <- OUT.LASSO
  
  PRED.ERROR.NORM      <- mean((Y1.test - X1.test%*%OUT.NORM$beta)^2)
  PRED.ERROR.COMB      <- mean((Y1.test - X1.test%*%OUT.COMB$beta)^2)
  PRED.ERROR.FIX       <- mean((Y1.test - X1.test%*%OUT.FIX$beta)^2)
  PRED.ERROR.LASSO     <- mean((Y1.test - X1.test%*%OUT.LASSO$beta)^2)
  
  OUTPUT$PRED.ERROR.true <- mean((Y1.test - X1.test%*%beta.x)^2)
  OUTPUT$PRED.ERROR.null <- mean((Y1.test)^2)
  
  
  return(OUTPUT = OUTPUT)
  
}

cv.par <- function(d, dt, crossval.ind, X1, Y1, M, beta.0, lambda.1, lambda.2, m.iter, n.iter, tol) {
  
  
  
  CVOUTPUT <-lapply(1:d, cv.loop,
                    dt, crossval.ind, X1, Y1, 
                    M, beta.0, 
                    lambda.1, lambda.2, m.iter, 
                    n.iter, tol)
  
  
  
  
  return(CVOUTPUT)
}

cv.loop <- function(ds, dt, crossval.ind, X1, Y1, M, beta.0, lambda.1, lambda.2, m.iter, n.iter, tol){
  
  CVOUTPUT <- NULL
  
  
  cat("Cross validation data set - ", ds, "\n")
  cat("Data set - ", dt, "\n")
  train.ind <- crossval.ind[[ds]]$train
  val.ind   <- crossval.ind[[ds]]$validate
  
  X1.train <- X1[train.ind,]
  Y1.train <- Y1[train.ind]
  
  X1.val   <- X1[val.ind,]
  Y1.val   <- Y1[val.ind]
  
  
  
  cat('*** Norm Lap ***')
  M.int               <- mat.to.laplacian(M, type = "normalized")
  output.train        <- lasso.net.grid(X1.train,Y1.train,beta.0,lambda.1,lambda.2, M.int, m.iter = m.iter, n.iter = n.iter, iscpp = T, tol = tol, alt.num = 4 )
  CVOUTPUT$beta.norm  <- output.train$beta
  
  cat('*** Fixed Norm ***')
  M.int               <- mat.to.laplacian(M, type = "normalized")  
  output.train        <- lasso.net.fixed(X1.train,Y1.train,beta.0,lambda.1,lambda.2, M.int, n.iter = n.iter, iscpp = T, tol = tol )
  CVOUTPUT$beta.fix   <- output.train$beta
  
  cat('*** Comb Lap ***')
  M.int               <- mat.to.laplacian(M, type = "combinatorial")
  output.train        <- lasso.net.grid(X1.train,Y1.train,beta.0,lambda.1,lambda.2, M.int, m.iter = m.iter, n.iter = n.iter, iscpp = T, tol = tol, alt.num = 4 )
  CVOUTPUT$beta.comb  <- output.train$beta
  
  cat('*** Lasso ***')
  output.train        <- lasso.net.grid(X1.train,Y1.train,beta.0,lambda.1,0, M, m.iter = m.iter, n.iter = n.iter, iscpp = T, tol = tol, alt.num = 4 )
  CVOUTPUT$beta.lasso <- output.train$beta
  
  
  CVOUTPUT$X1.val <- X1.val
  CVOUTPUT$Y1.val <- Y1.val 
  
  
  return(CVOUTPUT = CVOUTPUT)
}

CVfold <- function(n,K=10, is.random = FALSE){
  
  # determine if draw random indices or sequential
  if(is.random){
    id    <- sample(n,replace = FALSE)
  } else {
    id <- seq(n)  
  }
  
  # get indices:
  k <- as.integer(n*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),n),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  ind <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],validate=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(ind)
}











