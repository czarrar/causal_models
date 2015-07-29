#!/usr/bin/env Rscript

require(glasso)
require(cvTools)
require(lassoscore)

# note: should divide emp_cov by number of tests?
log_likelihood <- function(precision, emp_cov) {
  p      <- nrow(precision)
  logdet <- determinant(precision, logarithm=T)$modulus
  loglik <- 0.5 * (logdet - sum(emp_cov * precision) - p*log(2*pi))
  return(as.numeric(loglik))
}

# graphical lasso using cross-validation
glasso_cv <- function(ts, k=5, rholist=NULL, verbose=F) {
  if (is.null(rholist)) {
    # We will determine the maximum rho as one where the inverse covariance
    # has at least 5% non-zero off-diagonal elements (arbitrary)
    S       <- cov(ts)
    rholist <- seq(0, max(abs(S)), length=11)[-1]
    GLP     <- glassopath(S, rholist, trace=0)
    nonzeros<- apply(GLP$wi, 3, function(x) mean(x[upper.tri(x)]!=0))
    max.rho <- max(rholist[nonzeros > 0.05], rholist[1])
    
    # Now make the list of rhos
    rholist <- seq(0, max.rho, length=100)
  }
  
  n     <- nrow(ts)
  folds <- cvFolds(n, k, type="consecutive")
  
  loglikes <- sapply(1:k, function(ki) {
    if (verbose) cat("Fold ", ki, "\n")
    S_train <- cov(ts[folds$which!=ki,])
    S_test  <- cov(ts[folds$which==ki,])
    GLP     <- glassopath(S_train, rholist, trace=0)
    loglike <- apply(GLP$wi, 3, function(P_train) log_likelihood(P_train, S_test))
    loglike
  })
  
  ind     <- which.max(rowMeans(loglikes))
  rhomax  <- rholist[ind]
  S       <- cov(ts)
  a       <- glasso(S, rhomax)
  a$rhomax <- rhomax
  
  return(a)
}

#' Wrapper around Graphical Lasso
#'
#' Returns an adjacency matrix of connections. Does a cross-validation on the 
#' graphical lasso to determine the optimal lambda.
#'
#' @param ts_mat Time-series matrix (ntpts x nregions)
#' @param ... Other arguments/values to pass to the `glasso_cv` function
#'
#' @return matrix

conn_glassoscore <- function(ts_mat, scale=F, ret="p.sand", alpha=0.01, ...) {
  if (scale) ts_mat <- scale(ts_mat)
  tmp <- glasso_cv(ts_mat, ...)
  gls <- glassoscore(ts_mat, tmp$rhomax)
  adj <- (gls[[ret]] < alpha) * 1
  diag(adj) <- 0
  as.matrix(adj)
}

conn_mbscore <- function(ts_mat, scale=F, ret="p.sand", alpha=0.01, ...) {
  if (scale) ts_mat <- scale(ts_mat)
  tmp <- glasso_cv(ts_mat, ...)
  mbs <- mbscore(ts_mat, tmp$rhomax, family="gaussian")
  adj <- (mbs[[ret]] < alpha) * 1
  diag(adj) <- 0
  as.matrix(adj)
}

conn_cv_mbscore <- function(ts_mat, scale=F, ret="p.sand", alpha=0.01, ...) {
  if (scale) ts_mat <- scale(ts_mat)
  mbs <- cv.mbscore(ts_mat)
  adj <- (mbs[[ret]] < alpha) * 1
  diag(adj) <- 0
  as.matrix(adj)
}

# selects the appropriate lambda for the data
cv.mbscore <- function(x,k=10,subset=NULL,tol=1e-8,...){
	beta <- scores <- scorevar.mod <- scorevar.sand <- p.model <- p.sand <- matrix(NA,ncol(x),ncol(x))
	
	if(is.null(subset)) subset <- matrix(TRUE,ncol(x),ncol(x))
  
  n     <- nrow(x)
  folds <- cvFolds(n, k, type="consecutive")

	for(i in 1:ncol(x)){
    if (!(i %% 10)) cat(i,".",sep="")
    cv <- cv.glmnet(y=x[,i], x=x[,-i], foldid=folds$which, family="gaussian", 
                    thresh=tol, maxit=1000, standardize=FALSE)
		mod <- lassoscore(x[,i],x[,-i], lambda=cv$lambda.min, subset=which(subset[-i,i]), tol=tol, ...)
		scores[-i,i] <- mod$scores
		p.model[-i,i] <- mod$p.model
		p.sand[-i,i] <- mod$p.sand
		scorevar.mod[-i,i] <- mod$scorevar.model
		scorevar.sand[-i,i] <- mod$scorevar.sand
		beta[-i,i] <- mod$fit$beta	
	}
  cat("\n")
	re <- list("scores"=scores,
				"p.model" = p.model,
				"p.sand" = p.sand,
				"scorevar.mod" = scorevar.mod,
				"scorevar.sand" = scorevar.sand,
				"beta" = beta)
	re$lambda <- cv$lambda.min
	class(re) <- "mbscore"			
	return(re)
}

cv.mbscore2 <- function(x,k=10,subset=NULL,tol=1e-8,...){
	beta <- scores <- scorevar.mod <- scorevar.sand <- p.model <- p.sand <- matrix(NA,ncol(x),ncol(x))
	
	if(is.null(subset)) subset <- matrix(TRUE,ncol(x),ncol(x))
  
  n     <- nrow(x)
  folds <- cvFolds(n, k, type="consecutive")

	for(i in 1:ncol(x)){
    if (!(i %% 10)) cat(i,".",sep="")
    cv <- cv.glmnet(y=x[,i], x=x[,-i], foldid=folds$which, family="gaussian", 
                    thresh=tol, maxit=1000, standardize=FALSE, alpha=0.5)
		mod <- lassoscore2(x[,i],x[,-i], lambda=cv$lambda.min, subset=which(subset[-i,i]), tol=tol, ...)
		scores[-i,i] <- mod$scores
		p.model[-i,i] <- mod$p.model
		p.sand[-i,i] <- mod$p.sand
		scorevar.mod[-i,i] <- mod$scorevar.model
		scorevar.sand[-i,i] <- mod$scorevar.sand
		beta[-i,i] <- mod$fit$beta	
	}
  cat("\n")
	re <- list("scores"=scores,
				"p.model" = p.model,
				"p.sand" = p.sand,
				"scorevar.mod" = scorevar.mod,
				"scorevar.sand" = scorevar.sand,
				"beta" = beta)
	re$lambda <- cv$lambda.min
	class(re) <- "mbscore"			
	return(re)
}


#lassoscore <- function(y,X, lambda=0, family=c("gaussian","binomial","poisson"), tol = .Machine$double.eps, maxit=1000, 
#                       resvar = NULL, verbose=FALSE, subset = NULL){
#  family = match.arg(family)
#  family.obj <- get(family,mode="function",envir=parent.frame())()
#  
#  if(!(family %in% c("gaussian","binomial","poisson"))) stop("family not supported!")                   
#  X <- scale(X)*sqrt(nrow(X)/(nrow(X)-1))
#  if(family=="gaussian") y <- scale(y,scale=FALSE)
#  
#  ##initial fit
#  out0 <- glmnet(y=y,x=X, family=family, lambda=lambda, thresh = tol, maxit=maxit,standardize=FALSE)
#
#}


lassoscore2 <- function(y,X, lambda=0, family=c("gaussian","binomial","poisson"), tol = .Machine$double.eps, maxit=1000, 
                       resvar = NULL, verbose=FALSE, subset = NULL){
  family = match.arg(family)
  family.obj <- get(family,mode="function",envir=parent.frame())()
  
  if(!(family %in% c("gaussian","binomial","poisson"))) stop("family not supported!")                   
  X <- scale(X)*sqrt(nrow(X)/(nrow(X)-1))
  if(family=="gaussian") y <- scale(y,scale=FALSE)
  
  ##initial fit
  out0 <- glmnet(y=y,x=X, family=family, lambda=lambda, thresh = tol, maxit=maxit,standardize=FALSE,alpha=0.5)
  out0$r <- as.vector(family.obj$linkinv(as.vector(out0$a0 + X%*%out0$beta))-y)
  out0$v <- family.obj$var(family.obj$linkinv(as.vector(out0$a0 + X%*%out0$beta)))
  
  
  if(family =="gaussian" & !is.numeric(resvar)){
    resvar <- sum(out0$r^2)/(length(y)-sum(out0$beta !=0))
  } else if(family !="gaussian"){
    resvar <- 1
  }
  
  out0$n <- nrow(X)
  out0$beta <- as.vector(out0$beta)
  
  wh <- as.vector(out0$beta != 0)
  if(is.null(subset)){
    subset <- 1:ncol(X)
  }
  if(verbose){
    cat("\nProgress:\n")
    pb <- txtProgressBar(min = 0, max = length(subset), style = 3)
    pb.i <- 0
  }
  scores <- scorevar.sand.cons <- scorevar.model.cons <- scorevar.sand <- scorevar.model <- numeric(ncol(X))
  scores[-subset] <- NA
  
  for(i in subset){
  	if(out0$beta[i] != 0){
    	out <- glmnet(y=y,x=X[,-i], family=family, lambda=lambda, thresh = tol, maxit=maxit,standardize=FALSE,alpha=0.5)
    	out$r <- as.vector(family.obj$linkinv(out$a0 +as.vector(X[,-i]%*%out$beta))-y)
	    out$v <- family.obj$var(family.obj$linkinv(out$a0 +as.vector(X[,-i]%*%out$beta)))  		
    	out$n <- nrow(X)
    	out$beta <- as.vector(out$beta)
    	Xs <- X[,-i,drop=FALSE][,as.vector(out$beta !=0),drop=FALSE]
    } else {
    	out <- out0
    	Xs <- X[,as.vector(out$beta !=0),drop=FALSE]
    }
    
    xx <- X[,i,drop=FALSE]
    
    scores[i] <- sum(out$r*X[,i])/sqrt(out$n)
       
    scorevar.model.cons[i] <- with(out, resvar*(crossprod(xx*sqrt(v))/n))
    scorevar.sand.cons[i] <- with(out, var(xx*r)*(n-1)/n)

    if(ncol(Xs) == 0){
      scorevar.sand[i] <- scorevar.sand.cons[i]
      scorevar.model[i] <- scorevar.model.cons[i]	
    } else if(nrow(Xs) > ncol(Xs) & ncol(Xs) >0){ 
      Ui <- with(out, solve(crossprod(Xs*sqrt(v))/n))
      V <- with(out, var(r*Xs)*(n-1))
      Va <- with(out,cov(r*Xs,r*xx)*(n-1))
      Ua <-  with(out, crossprod(Xs,v*xx)/n)
      va <-  with(out, var(xx*r)*(n-1))
      
      scorevar.sand[i] <- (va + t(Ua)%*%(Ui)%*%(V%*%Ui%*%Ua - 2*Va))/with(out,n-sum(beta !=0))
      scorevar.model[i] <- with(out, resvar*(crossprod(xx*sqrt(v))/n - t(Ua)%*%Ui%*%Ua)	)  
    }
    if(verbose){
      pb.i <- pb.i+1
      setTxtProgressBar(pb, pb.i)
    }
  }
  re <- list(
    "fit"=out0,
    "scores" = scores,
    "scorevar.model.cons" = scorevar.model.cons,
    "scorevar.sand.cons" = scorevar.model.cons,
    "scorevar.model" = scorevar.model,
    "scorevar.sand" = scorevar.sand,	
    "p.model.cons" = pchisq(scores^2/scorevar.model.cons,df=1,lower.tail=FALSE),
    "p.sand.cons" = pchisq(scores^2/scorevar.sand.cons,df=1,lower.tail=FALSE),
    "p.model" = pchisq(scores^2/scorevar.model,df=1,lower.tail=FALSE),
    "p.sand" = pchisq(scores^2/scorevar.sand,df=1,lower.tail=FALSE),
    "lambda" = lambda)
  
  class(re) <- "lassoscore"
  if(verbose) close(pb)
  return(re)
}
