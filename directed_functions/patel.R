# See Ramsey et al., 2014 for psuedo-code
# and also Patel et al. ?

patel.tau <- function(x, S=NULL, to.scale=F, verbose=F) {
  if (to.scale) x <- scale(x)
  if (is.null(S)) S <- matrix(1, ncol(x), ncol(x))
  
  # 1. Map the data into the interval [0,1]:
  if (verbose) cat("map data into the interval [0,1]\n")
  X10 <- apply(x, 2, quantile, 0.1) # cutoff for the 10th percentile
  X90 <- apply(x, 2, quantile, 0.9) # cutoff for the 90th percentile
  a   <- sweep(sweep(x, 2, X10), 2, X90-X10, FUN="/") # (X-X10)/(X90-X10)
  x2  <- apply(a, c(1,2), function(v) max(min(v,1),0)) # keep in [0,1] interval
    
  # 2. theta1 = X2*Y2 or x2*t(x2)
  if (verbose) cat("theta1\n")
  theta1 <- crossprod(x2)

  # 3. theta2 = X2 * (1-Y2)
  if (verbose) cat("theta2\n")
  theta2 <- crossprod(x2, 1-x2)

  # 4. theta3 = (1-X2)*Y2
  if (verbose) cat("theta3\n")
  theta3 <- crossprod(1-x2,x2)

  # 5. Set tau
  if (verbose) cat("tau\n")
  tau  <- matrix(0, ncol(x), ncol(x))
  inds <- theta2 > theta3
  tau[inds] <- 1 - (theta1[inds] + theta3[inds])/(theta1[inds] + theta2[inds])
  tau[!inds] <- (theta1[!inds] + theta2[!inds])/(theta1[!inds] + theta3[!inds]) - 1

  # 6. Directed unweighted graph
  if (verbose) cat("dag\n")
  dag <- S * (tau<0)
    
  list(dag=dag, tau=tau)
}
