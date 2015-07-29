#!/usr/bin/env Rscript

library(huge)

## huge.npn(dat)
#system.time(tmp <- huge(dat, nlambda=30, method="ct", verbose=F))
#system.time(tmp <- huge(dat, nlambda=30, method="mb", verbose=F))
#system.time(tmp <- huge(dat, nlambda=30, method="glasso", verbose=F))
#
#res <- huge.select(tmp, "ebic")
#res <- huge.select(tmp, "stars")
#res <- huge.select(tmp, "ric")

# opts:
# expand.grid(list(npn=c(T,F), method=c("ct", "mb", "glasso"), scr=c(T,F), criterion=c("ric", "stars", "ebic")
# remove scr option T when method is ct
# remove criterion option ebic when method is ct or mb
# scr.num => n-1 or n/log(n)

conn_huge <- function(ts_mat, to.scale=F, 
                      npn=F, npn.func="shrinkage", 
                      method="ct", nlambda=40, scr=F, scr.num=NULL, 
                      criterion="ric", rep.num=30, 
                      verbose=F)
{
  if (to.scale) ts_mat <- scale(ts_mat)
  if (npn) ts_mat <- huge.npn(ts_mat, npn.func, verbose=verbose)
  
  # Estimate the solution path
  res.huge <- huge(ts_mat, nlambda=nlambda, method=method, scr=scr, verbose=verbose)
  
  # Select the graph
  res.select <- huge.select(res.huge, criterion=criterion, rep.num=rep.num, verbose=verbose)
  
  # For now only return the adjacencies
  as.matrix(res.select$refit)
}
