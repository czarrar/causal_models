#!/usr/bin/env Rscript

require(pcalg)
library(graph)

upper <- function(x) x[upper.tri(x)]

# Function to run Peter Clark Algorithm given (2D) time-series input
# the output is a 2D adjacency matrix

#' Wrapper around Peter Clark Algorithm
#'
#' Returns an adjacency matrix of connections. Calls the `skeleton` function 
#' in the `pcalg` package. The output is then converted to be a matrix.
#'
#' @param ts_mat Time-series matrix (ntpts x nregions)
#' @param indepTest Function for testing conditional independence (default: gaussCItest)
#' @param alpha Significance level for individual conditional independence tests (default: 0.01)
#' @param ... Other arguments/values to pass to the `skeleton` function
#'
#' @return matrix
conn_pc <- function(ts_mat, indepTest=gaussCItest, alpha=0.05, fdr.correct="none", ...) {
  suffStat  <- list(C = cor(ts_mat), n = nrow(ts_mat))
  sk.res    <- skeleton(suffStat, indepTest=indepTest, p=ncol(ts_mat), alpha=alpha, ...)
  
  if (fdr.correct == "qval") {
    pMat <- matrix(0, ncol(sk.res@pMax), ncol(sk.res@pMax))
    fvec <- fdrtool(upper(sk.res@pMax), statistic="pvalue", plot=F, verbose=F)
    pMat[upper.tri(pMat)] <- fvec$qval
    pMat <- pMat + t(pMat)
    diag(pMat) <- 1
  } else if (fdr.correct == "lfdr") {
    pMat <- matrix(0, ncol(sk.res@pMax), ncol(sk.res@pMax))
    fvec <- fdrtool(upper(sk.res@pMax), statistic="pvalue", plot=F, verbose=F)    
    pMat[upper.tri(pMat)] <- fvec$lfdr
    pMat <- pMat + t(pMat)    
    diag(pMat) <- 1
  } else if (fdr.correct == "fdr") {
    pMat <- matrix(0, ncol(sk.res@pMax), ncol(sk.res@pMax))
    fvec <- p.adjust(upper(sk.res@pMax), "fdr")
    pMat[upper.tri(pMat)] <- fvec
    pMat <- pMat + t(pMat)    
    diag(pMat) <- 1
  } else {
    pMat <- sk.res@pMax
  }
  
  adj_mat   <- (pMat < alpha)*1
  
  as.matrix(adj_mat)
}

conn_pc_pvals <- function(ts_mat, indepTest=gaussCItest, fdr.correct="none", ...) {
  suffStat  <- list(C = cor(ts_mat), n = nrow(ts_mat))
  sk.res    <- skeleton(suffStat, indepTest=indepTest, p=ncol(ts_mat), alpha=0.05, method="stable.fast", ...)
  
  if (fdr.correct == "qval") {
    pMat <- matrix(0, ncol(sk.res@pMax), ncol(sk.res@pMax))
    fvec <- fdrtool(upper(sk.res@pMax), statistic="pvalue", plot=F, verbose=F)
    pMat[upper.tri(pMat)] <- fvec$qval
    pMat <- pMat + t(pMat)
    diag(pMat) <- 1
  } else if (fdr.correct == "lfdr") {
    pMat <- matrix(0, ncol(sk.res@pMax), ncol(sk.res@pMax))
    fvec <- fdrtool(upper(sk.res@pMax), statistic="pvalue", plot=F, verbose=F)    
    pMat[upper.tri(pMat)] <- fvec$lfdr
    pMat <- pMat + t(pMat)    
    diag(pMat) <- 1
  } else if (fdr.correct == "fdr") {
    pMat <- matrix(0, ncol(sk.res@pMax), ncol(sk.res@pMax))
    fvec <- p.adjust(upper(sk.res@pMax), "fdr")
    pMat[upper.tri(pMat)] <- fvec
    pMat <- pMat + t(pMat)    
    diag(pMat) <- 1
  } else {
    pMat <- sk.res@pMax*0
    pMat[upper.tri(pMat)] <- upper(sk.res@pMax)
    pMat <- pMat + t(pMat)
    diag(pMat) <- 1
  }
  
  as.matrix(pMat)
}
