#!/usr/bin/env Rscript

require(glmnet)
require(cvTools)
require(bigmemory)

mbglmnet <- function(x, k=10, tol=1e-8, alpha=0.5, parallel=TRUE, verbose=FALSE, ...) {
  nr    <- nrow(x)
  nc    <- ncol(x)
  folds <- cvFolds(nr, k, type="consecutive")
  betas <- big.matrix(nc, nc, init=0, shared=TRUE)
  
  step <- floor(nc * 0.1) # 10% increments
  if (step == 0) step <- nc
	l_ply(1:nc, function(i) {
    if (!(i %% step) & verbose) cat(floor(i/step)/0.1,"% - ",sep="")
    cv <- suppressWarnings(
            cv.glmnet(y=x[,i], x=x[,-i], foldid=folds$which, family="gaussian", 
                        thresh=tol, maxit=1000, alpha=alpha, ...)
          )
    res <- glmnet(y=x[,i], x=x[,-i], family="gaussian", thresh=tol, 
                  lambda=cv$lambda.min, alpha=alpha, ...)
		betas[-i,i] <- as.numeric(res$beta)
  }, .parallel=parallel)
	if(verbose) cat("\n")
  
  #r     <- as.matrix(betas)
  #d     <- sqrt(abs(diag(betas)))
  #r     <- sweep(r, 1, d)
  #r     <- sweep(r, 2, d)
  
  as.matrix(betas)
}

conn_mbglmnet <- function(ts_mat, scale=F, alpha=NULL, ...) {
  betas <- mbglmnet(ts_mat, standardize=scale, ...)
  r   <- beta2cor(betas)
  adj <- (r != 0) * 1
  diag(adj) <- 0
  as.matrix(adj)
}

conn_mbglmnet_cors <- function(ts_mat, scale=F, alpha=NULL, ...) {
  betas <- mbglmnet(ts_mat, standardize=scale, ...)
  r     <- beta2cor(betas)
  r
}

conn_mbglmnet_scale <- function(ts_mat, scale=T, alpha=NULL, ...) {
  betas <- mbglmnet(ts_mat, standardize=scale, ...)
  r <- beta2cor(betas)
  adj <- (r != 0) * 1
  diag(adj) <- 0
  as.matrix(adj)
}

# borrowed from parcor (adalasso) package
beta2cor <- function(Beta, verbose=FALSE){
  Dummy <- Beta*t(Beta)
  if (verbose==TRUE){
    cat("\nNumber of pairwise regression coefficients with conflicting signs:", (sum((Dummy) < 0))/2, "\n")
    cat("Number of partial correlation coefficients greater than 1 in absolute value:",(sum((Dummy) >1))/2 ,
  "\n\n")
  }
  Dummy[Dummy<0] <- 0 # if a product is <0 the partial correlation coefficient is set to 0
  Dummy[Dummy>1] <- 1 # partial correlation coefficients should be in the range of [-1,1]
  P <- sign(Beta)*sqrt(Dummy)
  diag(P) <- 1
  return(P)
}
