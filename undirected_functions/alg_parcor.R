#!/usr/bin/env Rscript

library(parcor)

# Modified version of their function to do parallel computation and also output adalasso result
adalasso.net <- function(X, k=10,use.Gram=FALSE,verbose=FALSE,intercept=TRUE,parallel=TRUE)
{
  p <- ncol(X)
  
  X <- scale(X)     # data needs to be centered and standardized
  colnames(X) <- 1:p    # each column gets a name
                      
  B.adalasso <- matrix(0, nrow=p, ncol=p)
  
  if (verbose==TRUE) {
    cat(paste("Performing local (adaptive) lasso regressions\n"))
    cat(paste("Vertex no "))
  }
  
  # result: matrix with each row result of adalasso (row-wise)
  B.adalasso <- laply(1:p, function(i) {
    if (verbose==TRUE) {
      if (!(i %% 10)) cat(paste(i,"..."))
    }
    
    noti <- (1:p)[-i]
    yi <- X[ ,i]       ## response
    Xi <- X[ ,noti]    ## predicted by all other nodes with i missing
    res <- rep(0, p)
    
    ## perform adaptive lasso regression & extract regression coefficients  
    dummy <- adalasso(Xi, yi, k=k,use.Gram=use.Gram,both=TRUE,intercept=intercept)
    coefi.adalasso <- dummy$coefficients.adalasso
    res[-i] <- coefi.adalasso 
    
    res
  }, .parallel=parallel)
  colnames(B.adalasso) <- 1:p
  
  pcor.adalasso <- Beta2parcor(B.adalasso, verbose=verbose)
  
  cat(paste("\n"))
  
  return(pcor.adalasso)
}

conn_adalasso <- function(ts_mat, scale=F, thr=0, parallel=TRUE, ...) {
  if (scale) ts_mat <- scale(ts_mat)
  pcor <- adalasso.net(ts_mat, parallel=parallel, ...)
  adj <- (abs(pcor) > thr) * 1
  diag(adj) <- 0
  as.matrix(adj)
}
