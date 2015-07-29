# For each node, we look at all elements that are connected to that node. we ask for all combinations of those elements and regress them on our node, and ask what combination gives a residual with the most non-gaussian distribution
library(plyr)
library(bigmemory)

# inpired by: library(nortest)
r.NG <- function(x) 
{
  tol <- 1e-8
  x <- sort(x)
  n <- length(x)
  if (n < 8) stop("sample size must be greater than 7")
  p <- pnorm((x - mean(x))/sd(x))
  h <- (2 * seq(1:n) - 1) * (log(p + tol) + log(1 - rev(p) + tol))
  A <- -n - mean(h)
  return(A)
}

library(inline)
library(Rcpp)
# inpired by: library(nortest)
cat("...compiling")
cpp.NG <- cxxfunction(signature(x1="numeric"), 
                      plugin="Rcpp",
                      body='
NumericVector x2 = x1;
NumericVector x = clone(x2);
int n = x.size();
double tol = 1e-8;

std::sort(x.begin(), x.end());

NumericVector p = Rcpp::pnorm((x - mean(x))/sd(x));

NumericVector rev_p = clone(p);
std::reverse(rev_p.begin(), rev_p.end());
//::Rf_copyMostAttrib(p, rev_p); 

NumericVector h(n);
for (double i=0; i<n; i++) {
  h[i] = (2 * (i+1) - 1) * (log(p[i] + tol) + log(1 - rev_p[i] + tol));
}

double A = -n - mean(h);

return(Rcpp::wrap(A));
')
cat("...done\n")

#a <- rnorm(1000, mean=5, sd=3)
#all.equal(tester(a), cpp.NG(a))
#cat("\nNG\n")
#for (j in 1:10) print(system.time(for (i in 1:100) NG(a)))
#cat("\ncpp.NG\n")
#for (j in 1:10) print(system.time(for (i in 1:100) cpp.NG(a)))

#d1 <- rnorm(1000, mean = 5, sd = 3)
#d2 <- runif(1000, min = 2, max = 4)
#tester <- function(fun, ..., n=50) for (i in 1:n) fun(...)
#system.time(tester(ad.test, d1))
#system.time(tester(ad.test, d2))
#system.time(tester(NG, d1))
#system.time(tester(NG, d2))
#NG <- function(x) as.numeric(ad.test(x)$statistic)

calc_residuals <- function(y, X) {
  if (!is.matrix(y))
    y <- as.matrix(y)
  if (!is.matrix(X))
    X <- as.matrix(X)
  
  # beta values
  iXtX <- solve(t(X) %*% X)
  b <- iXtX %*% t(X) %*% y
  
  # residuals
  res <- y - X %*% b
  
  res
}

calc_residuals_cors <- function(y, X, rs=NULL) {
  if (!is.matrix(y))
    y <- as.matrix(y)
  if (!is.matrix(X))
    X <- as.matrix(X)
  
  # betas / correlations
  if (is.null(rs)) rs <- t(cor(y, X))
  
  # residuals
  res <- y - X %*% rs
  
  res
}



#' Rule 1 for LOFS
#'
#' TODO
#'
#' @param S adjacency matrix of all connections between nodes (nnodes x nnodes)
#' @param dat matrix of time-series (ntpts x nnodes)
#' @return G matrix with connections in S oriented
#'
#' @references
#' TODO
#'
lofs.r1 <- function(S, dat, .NG=cpp.NG) {  
  nnodes  <- ncol(dat)
  
  # standardize
  dat     <- scale(dat)
  
  # output matrix
  G       <- matrix(0, nnodes, nnodes)
  
  # no self-nodes
  diag(S) <- 0
  
  for (i in 1:nnodes) {
    possible_parents <- which(S[,i]>0)
    
    # want to get every possible combination of parent nodes
    combos <- llply(1:length(possible_parents), function(j) {
      combn(possible_parents, j, simplify=F)
    })
    combos <- unlist(combos, recursive=F)
    
    # regress on node with possible parent combos
    # and get the non-gaussianity of the residuals
    ngs <- sapply(combos, function(inds) {
      resids <- calc_residuals(dat[,i,drop=F], dat[,inds,drop=F])
      .NG(resids)
    })
    
    # get non-gaussianity without anything
    self_ng <- .NG(dat[,i])
    
    # which node has the maximum non-gaussianity?
    # and greater than empty set
    max_ind <- which.max(ngs)
    if (ngs[max_ind] > self_ng) {
      G[combos[[max_ind]],i] <- S[combos[[max_ind]],i]
    }
  }
  
  return(G)
}

#' Rule 2 for LOFS
#' 
#' Basically for this one, I want to take the combo approach from lofs1
#' except that this time would want to do both directions X<-Y and Y->X
#' save each of those scores.... calculation for the score shoud be looked up
lofs.r2 <- function(S, dat, .NG=cpp.NG) {
  nnodes  <- ncol(dat)
  
  # standardize
  dat     <- scale(dat)
  
  # no self-nodes in input adjacency
  diag(S) <- 0
  
  
}

#' Rule 3 for LOFS
#' note X->Y is when M[i,j] == 1 (lower to higher numbered nodes)
lofs.r3 <- function(S, dat, rmat=NULL, to.scale=T, .NG=cpp.NG, 
                    fast=T, verbose=F, parallel=T) {
  if (verbose) {
    progress <- "text"
  } else {
    progress <- "none"
  }
  nnodes  <- ncol(dat)
  
  # standardize
  if (to.scale) dat <- scale(dat)
    
  # get correlation for initializing parameter search
  if (is.null(rmat)) rmat <- cor(dat)
  
  # no self-nodes in input adjacency
  diag(S) <- 0
  
  if (fast) {
    # NGs of individual time-series
    # NG(X) or NG(Y)
    node_NGs <- aaply(dat, 2, .NG, .parallel=parallel, .progress=progress)
  
    # Y = BX + E and get NG(E)
    # X <- Y
    # if NG(X,Y) + NG(Y) > NG(Y,X) + NG(X)
    # X -> Y
    # if NG(X,Y) + NG(Y) < NG(Y,X) + NG(X)
    sNG_mat <- big.matrix(nnodes, nnodes, init=0, shared=TRUE)
    l_ply(1:nnodes, function(ni) {
      adj_ni <- which(S[ni,]==1)
      resids <- calc_residuals_cors(dat[,adj_ni], dat[,ni], rmat[ni,adj_ni])  # NG(Y,X)
      toY    <- apply(resids, 2, .NG) + node_NGs[ni]    # NG(Y,X) + NG(X)
      sNG_mat[ni,adj_ni] <- toY
    }, .parallel=parallel, .progress=progress)
    sNG_mat <- as.matrix(sNG_mat)
  
    # If do a transpose, can then do comparison
    G <- (sNG_mat > t(sNG_mat)) * 1
  } else {
    # output matrix
    G       <- big.matrix(nnodes, nnodes, init=0, shared=T)
    
    # get i,j coordinates for edges
    inds    <- which(upper.tri(S) & S==1)
    coords  <- expand.grid(list(i=1:nnodes, j=1:nnodes))
    coords  <- coords[inds,]
    
    # Loop through each edge
    #for (ii in seq_along(inds)) {
    l_ply(seq_along(inds), function(ii) {
      i <- coords$i[ii]; j <- coords$j[ii]
      X <- dat[,i]; Y <- dat[,j]
      r <- rmat[i,j]
      
      toX <- .NG(calc_residuals_cors(X,Y,r)) + .NG(Y)
      toY <- .NG(calc_residuals_cors(Y,X,r)) + .NG(X)
      if (toX > toY) {
        G[i,j] <- 1
      } else {
        G[j,i] <- 1
      }
    }, .parallel=parallel, .progress=progress)
    
    G <- t(as.matrix(G))
  }
  
  G
}

# i think connections are oriented j->i so x[j,i] that's pos is j -> i
# but at the end, note X->Y is when M[i,j] == 1 (lower to higher numbered nodes)
#' Rule 4 for LOFS
#'
#' TODO
#'
#' @param S adjacency matrix of all connections between nodes (nnodes x nnodes)
#' @param dat matrix of time-series (ntpts x nnodes)
#' @param epsilon threshold for weighted matrix (default=0.1)
#' @param zeta range for free parameter, -zeta to zeta (default=1)
#' @param .NG function that measures the non-gaussianity of the data (default=NG aka the anderson-darling test)
#' @param add.self whether to have 1s in the diagonal of the weighted matrix (default=TRUE)
#' @return list(W,G) W weighted and G unweighted matrix with connections in S oriented
#'
#' @references
#' TODO
#'
lofs.r4 <- function(S, xdat, to.scale=T, epsilon=0.2, zeta=2, cordat=NULL, 
                    .NG=cpp.NG, verbose=FALSE, parallel=FALSE) 
{
  if (verbose) {
    suppressMessages(library(niftir)) # for progressbar
    vcat <- function(msg, ...) cat(sprintf(msg, ...), "\n")
  } else {
    vcat <- function(x) invisible(NULL)
  }
  
  nnodes  <- ncol(xdat)
  
  # standardize
  if (to.scale) xdat <- scale(xdat)
  
  # get correlation for initializing parameter search
  if (is.null(cordat)) cordat <- crossprod(xdat)/(nrow(xdat)-1)
  
  # no self-nodes in input adjacency
  diag(S) <- 0
  
  # construct W - output matrix
  W       <- matrix(0, nnodes, nnodes)
  diag(W) <- 1
  W[S!=0] <- NA
  ## save for now as big matrix
  W <- as.big.matrix(W, shared=TRUE)
  
  # function to be applied to each row
  # @param Wi = free parameters ? weights
  # @param X = self node time-series and adjacent time-series
  fun <- function(Wi, X) {
    .NG(tcrossprod(c(1,Wi), X))
  }
  
  select_nodes <- which(rowSums(S) > 0)
  select_nnodes <- length(select_nodes)
  if (verbose) cat("selecting ", select_nnodes, "/", nnodes, " nodes\n", sep="")
  
  # Loop through each row  
  if (verbose) pb <- progressbar(select_nnodes)
  #for (i in 1:nnodes) {
  l_ply(1:select_nnodes, function(ii) {
    i         <- select_nodes[ii]
    inds      <- as.integer(which(S[i,] > 0))
    # determine weights for free parameters while maximizing NG
    init.vals <- as.numeric(cordat[i,inds])
    #init.vals <- rep(0, length(inds))
    X         <- xdat[,c(i,inds),drop=F]
    res       <- optim(init.vals, fun, X=X, 
                       method="L-BFGS-B", control=list(fnscale=-1), 
                       lower=-zeta, upper=zeta)
    # assign optimized weights
    W[i,inds] <- res$par
    
    if (verbose && !(ii %% 10)) update(pb, ii)
  }, .parallel=parallel)
  if (verbose) end(pb)
  
  # Invert the matrix (TODO: make cleaner)
  W <- as.matrix(W)
  W <- t(W)
  
  # Threshold W to get G
  G <- (abs(W) > epsilon) * 1 # not sure if I use the absolute value
  
  # Other way to orient
  dag <- (W > t(W))*1
  
  list(W=W, G=G, dag=dag)
}

