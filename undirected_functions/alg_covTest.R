library(covTest)

conn_covtest <- function(ts_mat, alpha=0.01) {
  p.vals <- matrix(NA, ncol(ts_mat), ncol(ts_mat))
  for (i in 1:ncol(ts_mat)) {
    x <- ts_mat[,-i]; y <- ts_mat[,i]
    a <- lars(x,y)
    p.vals[-i,i] <- covTest(a, x, y)$results[,3]
  }
  adj <- (p.vals < alpha) * 1
  diag(adj) <- 0
  as.matrix(adj)
}


conn_covtest2 <- function(ts_mat, alpha=0.05, lambda2=0.5) {
  p.vals <- matrix(NA, ncol(ts_mat), ncol(ts_mat))
  for (i in 1:ncol(ts_mat)) {
    x <- ts_mat[,-i]; y <- ts_mat[,i]
    a <- lars.en(x,y,lambda2=lambda2)
    p.vals[-i,i] <- covTest(a, x, y)$results[,3]
  }
  adj <- (p.vals < alpha) * 1
  diag(adj) <- 0
  as.matrix(adj)
}


# other stuff

mbscorev2 <- function (x, subset = NULL, tol = 1e-08, ...) 
{
  beta <- scores <- scorevar.mod <- scorevar.sand <- p.model <- p.sand <- matrix(NA, 
                                                                                 ncol(x), ncol(x))
  if (is.null(subset)) 
    subset <- matrix(TRUE, ncol(x), ncol(x))
  for (i in 1:ncol(x)) {
    cvg <- cv.glmnet(x[,-i], x[,i])
    lambda <- cvg$lambda.min
    mod <- lassoscore(x[, i], x[, -i], lambda = lambda, subset = which(subset[-i, 
                                                                              i]), tol = tol, ...)
    scores[-i, i] <- mod$scores
    p.model[-i, i] <- mod$p.model
    p.sand[-i, i] <- mod$p.sand
    scorevar.mod[-i, i] <- mod$scorevar.model
    scorevar.sand[-i, i] <- mod$scorevar.sand
    beta[-i, i] <- mod$fit$beta
  }
  re <- list(scores = scores, p.model = p.model, p.sand = p.sand, 
             scorevar.mod = scorevar.mod, scorevar.sand = scorevar.sand, 
             beta = beta)
  re$lambda <- lambda
  class(re) <- "mbscore"
  return(re)
}

conn_mbscorev2 <- function(ts_mat, scale=F, ret="p.sand", alpha=0.01) {
  if (scale) ts_mat <- scale(ts_mat)
  mbs <- mbscorev2(ts_mat, family="gaussian")
  adj <- (mbs[[ret]] < alpha) * 1
  diag(adj) <- 0
  as.matrix(adj)
}
