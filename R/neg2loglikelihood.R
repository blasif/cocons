
#' GetNeg2loglikelihoodTaper
#' @description compute the negative 2 log likelihood based on theta
#' 
#' @usage GetNeg2loglikelihoodTaper(theta, par.pos, ref_taper, locs, 
#' x_covariates, smooth.limits, cholS, z, n, lambda, safe = TRUE)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list from \link{getDesignMatrix}.
#' @param ref_taper \code{(S4)} spam object based on a compact-supported covariance function.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param cholS \code{(S4)} Cholesky object from spam.
#' @param z \code{(numeric vector)} a vector of observed values.
#' @param n \code{(numeric)} dim(z)\[1\].
#' @param lambda \code{(numeric vector)} a vector including lambda.Sigma, lambda.betas, lambda.reg (in this order)
#' @param safe \code{(TRUE/FALSE)} if \code{TRUE} returns a large pre-defined value under Cholesky errors. Default \code{TRUE}.
#' @returns value 
#' @author Federico Blasi
GetNeg2loglikelihoodTaper <- function(theta, par.pos, ref_taper, locs,
                                      x_covariates, smooth.limits,
                                      cholS, z, n, lambda, safe = TRUE) {
  
  theta_list <- getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  ref_taper@entries <- ref_taper@entries * cov_rns_taper(theta = theta_list[-1], 
                                                                         locs = locs, 
                                                                         x_covariates =  x_covariates, 
                                                                         colindices = ref_taper@colindices, 
                                                                         rowpointers = ref_taper@rowpointers,
                                                                         smooth_limits =  smooth.limits)
  
  cholS <- tryCatch(spam::update.spam.chol.NgPeyton(cholS, ref_taper), error = function(e) e)

  if (inherits(cholS, "error")) {
    if(safe){return(1e+06)} else{
      stop("Cholesky error")
    }
  }
  
  logdet <- c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus)
  
  sumlogs <- sum(apply(z, 2, function(x){
    
    resid <- x - c(x_covariates %*% theta_list$mean)
    
    n * log(2 * pi) + 2 * c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) + 
      c(crossprod(spam::forwardsolve.spam(cholS,resid)))
  }))
  
  return(sumlogs + .cocons.getPen(n * dim(z)[2], lambda, theta_list, smooth.limits))
  
}

#' GetNeg2loglikelihoodTaperProfile
#' @description compute the negative 2 log likelihood based on theta
#' 
#' @usage GetNeg2loglikelihoodTaperProfile(theta, par.pos, ref_taper, 
#' locs, x_covariates, smooth.limits, cholS, z, n, lambda, safe = TRUE)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list.
#' @param ref_taper \code{(S4)} spam object based on a taper based covariance function.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param cholS \code{(S4)} Cholesky object from spam.
#' @param z \code{(numeric vector)} a vector of observed values.
#' @param n \code{(integer)} dim(z)\[1\].
#' @param lambda \code{(numeric vector)} a vector including lambda.Sigma, lambda.betas, lambda.reg (in this order)
#' @param safe \code{(TRUE/FALSE)} if \code{TRUE} returns a large pre-defined value under Cholesky errors. Default \code{TRUE}.
#' @returns \code{(numeric)} 
#' @author Federico Blasi
GetNeg2loglikelihoodTaperProfile <- function(theta, par.pos, ref_taper, locs,
                                             x_covariates, smooth.limits,
                                             cholS, z, n, lambda, safe = TRUE) {
  
  theta_list <- getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  theta_list$std.dev[1] <- 0
  
  ref_taper@entries <- ref_taper@entries * cov_rns_taper(theta = theta_list[-1], 
                                                                         locs = locs, 
                                                                         x_covariates =  x_covariates, 
                                                                         colindices = ref_taper@colindices, 
                                                                         rowpointers = ref_taper@rowpointers,
                                                                         smooth_limits =  smooth.limits)
  
  cholS <- tryCatch(spam::update.spam.chol.NgPeyton(cholS, ref_taper), error = function(e) e)
  
  if (inherits(cholS, "error")) {
    if(safe){return(1e+06)} else{
      stop("Cholesky error")
    }
  }
  
  logdet <- c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus)
  
  sum_in <- sum(apply(z,2,function(x){
    c(crossprod(spam::forwardsolve.spam(cholS, x - c(x_covariates %*% theta_list$mean))))
  }))
  
  return(dim(z)[2] * n * log(2 * pi) + 
           dim(z)[2] * n + 
           dim(z)[2] * 2 * logdet + 
           dim(z)[2] * n * log(sum_in / (dim(z)[2] * n)) + 
           .cocons.getPen(n * dim(z)[2], lambda, theta_list, smooth.limits))

}

#' GetNeg2loglikelihoodProfile
#' @description compute the negative 2 log likelihood based on theta
#'
#' @usage GetNeg2loglikelihoodProfile(theta, par.pos, locs, x_covariates, 
#' smooth.limits, z, n, x_betas,lambda, safe = TRUE)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param z \code{(numeric vector)} a vector of observed values.
#' @param n \code{(integer)} dim(z)\[1\].
#' @param x_betas \code{(matrix) or (data.frame)} design matrix for the spatial mean.
#' @param lambda \code{(numeric vector)} a vector including lambda.Sigma, lambda.betas, lambda.reg (in this order)
#' @param safe \code{(TRUE/FALSE)} if \code{TRUE} returns a large pre-defined value under Cholesky errors. Default \code{TRUE}.
#' @returns value 
#' @author Federico Blasi
GetNeg2loglikelihoodProfile <- function(theta, par.pos, locs, x_covariates, 
                                        smooth.limits, z, n, x_betas, lambda, safe = TRUE) {
  
  theta_list <- cocons::getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  Sigma_cpp <- cocons::cov_rns(theta = theta_list[-1], locs = locs, 
                              x_covariates =  x_covariates,
                              smooth_limits = smooth.limits)
  
  cholS <- tryCatch(base::chol(Sigma_cpp), error = function(e) e)
  
  if (inherits(cholS, "error")) {
    if(safe){return(1e+06)} else{
      stop("Cholesky error")
    }
  }
  
  if(T){
    V <- backsolve(cholS, forwardsolve(cholS, x_betas,
                                       transpose = TRUE, 
                                       upper.tri = TRUE))
    W <- crossprod(x_betas, V)
    Sigma_inv <- chol2inv(cholS)
    P_mat <- Sigma_inv - V %*% solve(W, t(V))
  }
  
  logdet <- sum(log(diag(cholS)))
  
  sum_logliks <- sum(apply(z, 2, function(x){
    
    quad_sum <- c(crossprod(x, P_mat %*% x))
    
    return(n * log(2 * pi) + 2 * logdet + quad_sum)
  }))
  
  return(sum_logliks + 
           .cocons.getPen(n * dim(z)[2], lambda, theta_list, smooth.limits))
  
}
    
#' GetNeg2loglikelihood
#' @description compute the negative 2 log likelihood based on theta
#'
#' @usage GetNeg2loglikelihood(theta, par.pos, locs, x_covariates, 
#' smooth.limits, z, n, lambda, safe = TRUE)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param z \code{(numeric vector)} a vector of observed values.
#' @param n \code{(integer)} dim(z)\[1\].
#' @param lambda \code{(numeric vector)} a vector including lambda.Sigma, lambda.betas, lambda.reg (in this order)
#' @param safe \code{(TRUE/FALSE)} if \code{TRUE} returns a large pre-defined value under Cholesky errors. Default \code{TRUE}.
#' @returns value
#' @author Federico Blasi
GetNeg2loglikelihood <- function(theta, 
                                 par.pos, 
                                 locs, 
                                 x_covariates, 
                                 smooth.limits, 
                                 z, 
                                 n,
                                 lambda,
                                 safe = TRUE) {
  
  theta_list <- cocons::getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  Sigma_cpp <- cocons::cov_rns(theta = theta_list[-1], 
                               locs = locs, 
                               x_covariates =  x_covariates,
                               smooth_limits = smooth.limits)
  
  cholS <- tryCatch(base::chol(Sigma_cpp), error = function(e) e)
  
  if (inherits(cholS, "error")) {
    if(safe){return(1e+06)} else{
      stop("Cholesky error")
    }
  }
  
  logdet <- sum(log(diag(cholS)))
  sum_logliks <- 0
  tmp_trend <- c(x_covariates %*% theta_list$mean)
  
  sum_logliks <- sum(apply(z, 2, function(x){
    tmp_a <- x - tmp_trend
    n * log(2 * pi) + 2 * logdet + c(crossprod(forwardsolve(cholS, 
                                                        tmp_a, 
                                                        transpose = TRUE, 
                                                        upper.tri = TRUE)))
  }))
  
  return(sum_logliks +  .cocons.getPen(n * dim(z)[2], lambda, theta_list, smooth.limits))
  
}

#' GetNeg2loglikelihoodREML
#' @description compute the negative 2 log REML likelihood based on theta
#'
#' @usage GetNeg2loglikelihoodREML(theta, par.pos, locs, x_covariates, x_betas,
#' smooth.limits, z, n, lambda, safe = TRUE)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param x_betas \code{(matrix) or (data.frame)} design matrix for the spatial mean.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param z \code{(numeric vector)} a vector of contrasts.
#' @param n \code{(integer)} dim(z)\[1\].
#' @param lambda \code{(numeric vector)} a vector including lambda.Sigma, lambda.betas, lambda.reg (in this order)
#' @param safe \code{(TRUE/FALSE)} if \code{TRUE} returns a large pre-defined value under Cholesky errors. Default \code{TRUE}.
#' @returns value
#' @author Federico Blasi
GetNeg2loglikelihoodREML <- function(theta, 
                                      par.pos, 
                                      locs, 
                                      x_covariates, 
                                      x_betas,
                                      smooth.limits, 
                                      z, 
                                      n,
                                      lambda,
                                      safe = TRUE) {
  
  theta_list <- cocons::getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  Sigma_cpp <- cocons::cov_rns(theta = theta_list[-1], 
                               locs = locs, 
                               x_covariates =  x_covariates,
                               smooth_limits = smooth.limits)
  
  cholS <- tryCatch(base::chol(Sigma_cpp), error = function(e) e)
  
  if (inherits(cholS, "error")) {
    if(safe){return(1e+06)} else{
      stop("Cholesky error")
    }
  }
  
  logdet <- sum(log(diag(cholS)))
  sum_logliks <- 0

  p <- qr(x_covariates)$rank
  
  if(T){
    V <- backsolve(cholS, forwardsolve(cholS, x_covariates,
                                       transpose = TRUE, 
                                       upper.tri = TRUE))
    W <- crossprod(x_covariates, V)
    Sigma_inv <- chol2inv(cholS)
    P_mat <- Sigma_inv - V %*% solve(W, t(V))
    
    chol_W <- chol(W)
  }
  
  sum_logliks <- sum(apply(z, 2, function(x){
    
    (n - p) * log(2 * pi) + 2 * logdet + 2 * sum(log(diag(chol_W))) +  c(crossprod(x, P_mat %*% x))
    
  }))
  
  return(sum_logliks +  .cocons.getPen((n - p) * dim(z)[2], lambda, theta_list, smooth.limits))
  
}



