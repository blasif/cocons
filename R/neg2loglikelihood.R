
#' GetNeg2loglikelihoodTaper
#' @description compute the negative 2 log likelihood based on theta
#' 
#' @usage GetNeg2loglikelihoodTaper(theta, par.pos, ref_taper, locs, 
#' x_covariates, smooth.limits, cholS, z, n, lambda)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list from \link{getDesignMatrix}.
#' @param ref_taper \code{(S4)} spam object based on a compact-supported covariance function.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param cholS \code{(S4)} Cholesky object from spam.
#' @param z \code{(numeric vector)} a vector of observed values.
#' @param n \code{(numeric)} dim(z)\[1\].
#' @param lambda \code{(numeric)} regularization parameter.
#' @returns value 
#' @author Federico Blasi
GetNeg2loglikelihoodTaper <- function(theta, par.pos, ref_taper, locs,
                                      x_covariates, smooth.limits,
                                      cholS, z, n, lambda) {
  
  theta_list <- getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  ref_taper@entries <- ref_taper@entries * cov_rns_taper(theta = theta_list[-1], 
                                                                         locs = locs, 
                                                                         x_covariates =  x_covariates, 
                                                                         colindices = ref_taper@colindices, 
                                                                         rowpointers = ref_taper@rowpointers,
                                                                         smooth_limits =  smooth.limits)
  
  Sigma_cpp <- spam::update.spam.chol.NgPeyton(cholS, ref_taper)
  logdet <- c(spam::determinant.spam.chol.NgPeyton(Sigma_cpp)$modulus)
  
  sumlogs <- 0
  
  for(ii in 1:dim(z)[2]){
    
    resid <- z[,ii] - c(x_covariates %*% theta_list$mean)
    
    sumlogs <- sumlogs + n * log(2 * pi) + 2 * c(spam::determinant.spam.chol.NgPeyton(Sigma_cpp)$modulus) + 
      sum(resid * spam::solve.spam(Sigma_cpp, resid))
    
  }
  
  return(sumlogs + getPen(n * dim(z)[2], lambda, theta_list, smooth.limits))

}

#' GetNeg2loglikelihoodTaperProfile
#' @description compute the negative 2 log likelihood based on theta
#' 
#' @usage GetNeg2loglikelihoodTaperProfile(theta, par.pos, ref_taper, 
#' locs, x_covariates, smooth.limits, cholS, z, n, lambda)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list.
#' @param ref_taper \code{(S4)} spam object based on a taper based covariance function.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param cholS \code{(S4)} Cholesky object from spam.
#' @param z \code{(numeric vector)} a vector of observed values.
#' @param n \code{(integer)} dim(z)\[1\].
#' @param lambda \code{(numeric)} regularization parameter.
#' @returns value 
#' @author Federico Blasi
GetNeg2loglikelihoodTaperProfile <- function(theta, par.pos, ref_taper, locs,
                                             x_covariates, smooth.limits,
                                             cholS, z, n, lambda) {
  
  theta_list <- getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  theta_list$std.dev[1] <- 0
  
  ref_taper@entries <- ref_taper@entries * cov_rns_taper(theta = theta_list[-1], 
                                                                         locs = locs, 
                                                                         x_covariates =  x_covariates, 
                                                                         colindices = ref_taper@colindices, 
                                                                         rowpointers = ref_taper@rowpointers,
                                                                         smooth_limits =  smooth.limits)
  
  Sigma_cpp <- spam::update.spam.chol.NgPeyton(cholS, ref_taper)
  logdet <- c(spam::determinant.spam.chol.NgPeyton(Sigma_cpp)$modulus)
  
  resid <- z - c(x_covariates %*% theta_list$mean)
  
  return(n * log(2 * pi / n) + n + 2 * logdet + 
    n * log(sum(resid * spam::solve.spam(Sigma_cpp, resid))) + getPen(n, lambda, theta_list, smooth.limits))

}

#' GetNeg2loglikelihoodProfile
#' @description compute the negative 2 log likelihood based on theta
#'
#' @usage GetNeg2loglikelihoodProfile(theta, par.pos, locs, x_covariates, 
#' smooth.limits, z, n, x_betas,lambda)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param z \code{(numeric vector)} a vector of observed values.
#' @param n \code{(integer)} dim(z)\[1\].
#' @param x_betas \code{(matrix) or (data.frame)} design matrix for the trend.
#' @param lambda \code{(numeric)} regularization parameter.
#' @returns value 
#' @author Federico Blasi
GetNeg2loglikelihoodProfile <- function(theta, par.pos, locs, x_covariates, 
                                        smooth.limits, z, n, x_betas, lambda) {
  
  theta_list <- cocons::getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  Sigma_cpp <- cocons::cov_rns(theta = theta_list[-1], locs = locs, 
                              x_covariates =  x_covariates,
                              smooth_limits = smooth.limits)
  
  # Compute P
  SigmaX <- solve(Sigma_cpp, x_betas)
  P_tete <- solve(Sigma_cpp) - SigmaX %*% solve(t(x_betas) %*% SigmaX) %*% t(SigmaX)
  
  check_pd <- tryCatch(cholS <- base::chol(Sigma_cpp), error = function(e) e)
  if (inherits(check_pd, "error")) {
    return(1e+06)
  } else{
    
    logdet <- sum(log(diag(cholS)))
    quad_sum <- t(z) %*% P_tete %*% z
    
    temp_thetas <- cocons::getModelLists(theta = theta, par.pos = par.pos, type = "classic")
    
    return(n * log(2 * pi) + 2 * logdet + quad_sum + 
      getPen(n, lambda, theta_list, smooth.limits))
    
  }
}
    
#' GetNeg2loglikelihood
#' @description compute the negative 2 log likelihood based on theta
#'
#' @usage GetNeg2loglikelihood(theta, par.pos, locs, x_covariates, 
#' smooth.limits, z, n, lambda)
#' @param theta \code{(numeric vector)} a vector with parameters values.
#' @param par.pos \code{(list)} par.pos list.
#' @param locs \code{(matrix)} spatial location matrix.
#' @param x_covariates \code{(data.frame)} design matrix.
#' @param smooth.limits \code{(numeric vector)} smooth.limits.
#' @param z \code{(numeric vector)} a vector of observed values.
#' @param n \code{(integer)} dim(z)\[1\].
#' @param lambda \code{(numeric)} regularization parameter.
#' @returns value
#' @author Federico Blasi
GetNeg2loglikelihood <- function(theta, 
                                 par.pos, 
                                 locs, 
                                 x_covariates, 
                                 smooth.limits, 
                                 z, 
                                 n,
                                 lambda) {
  
  theta_list <- cocons::getModelLists(theta = theta, par.pos = par.pos, type = "diff")
  
  Sigma_cpp <- cocons::cov_rns(theta = theta_list[-1], 
                               locs = locs, 
                               x_covariates =  x_covariates,
                               smooth_limits = smooth.limits)
  
  check_pd <- tryCatch(cholS <- base::chol(Sigma_cpp), error = function(e) e)
  if (inherits(check_pd, "error")) {
    return(1e+06)
  } else{
    
    logdet <- sum(log(diag(cholS)))
    sum_logliks <- 0
    
    for(ii in 1:dim(z)[2]){
      
      resid <- c(z[,ii]) - c(x_covariates %*% theta_list$mean)
      
      sum_logliks <- sum_logliks + n * log(2 * pi) + 2 * logdet + sum(resid * backsolve(cholS,
                                                                                        forwardsolve(cholS, 
                                                                                                     resid, 
                                                                                                     transpose = TRUE, 
                                                                                                     upper.tri = TRUE),
                                                                                        n)) 
      
    }
    
    return(sum_logliks +  getPen(n * dim(z)[2], lambda, theta_list, smooth.limits))
    
  }
}
