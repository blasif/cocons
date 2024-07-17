
#' GetNeg2loglikelihoodTaper
#' @description compute the negative 2 log likelihood based on theta
#' 
#' @usage function(theta)
#' @param theta a vector with parameters values
#' @param par.pos a vector with parameters values
#' @param ref_taper a vector with parameters values
#' @param locs a vector with parameters values
#' @param x_covariates a vector with parameters values
#' @param smooth.limits a vector with parameters values
#' @param cholS a vector with parameters values
#' @param z a vector with parameters values
#' @param n a vector with parameters values
#' @returns neg2loglikelihood 
#' @examples 
#' rnorm(1)
#' @author Federico Blasi
GetNeg2loglikelihoodTaper <- function(theta, par.pos, ref_taper, locs,
                                      x_covariates, smooth.limits,
                                      cholS, z, n, lambda) {
  
  theta_list <- getModelLists(theta = theta, par.pos = par.pos, type = 'diff')
  
  ref_taper@entries <- ref_taper@entries * cov_rns_taper_optimized_range(theta = theta_list[-1], 
                                                                         locs = locs, 
                                                                         x_covariates =  x_covariates, 
                                                                         colindices = ref_taper@colindices, 
                                                                         rowpointers = ref_taper@rowpointers,
                                                                         smooth_limits =  smooth.limits)
  
  Sigma_cpp <- spam::update.spam.chol.NgPeyton(cholS, ref_taper)
  logdet <- c(determinant.spam.chol.NgPeyton(Sigma_cpp)$modulus)
  
  resid <- z - c(x_covariates %*% theta_list$mean)
  
  return(n * log(2 * pi) + 2 * c(determinant.spam.chol.NgPeyton(Sigma_cpp)$modulus) + 
           sum(resid * solve.spam(Sigma_cpp, resid)) + 
           n * 2 * lambda * exp(theta_list$scale[1]) * 
           sqrt(((smooth.limits[2]-smooth.limits[1])/ 
                       (1 + exp(-theta_list$smooth[1])) + 
                  smooth.limits[1]))
  )
}

GetNeg2loglikelihoodTaperProfile <- function(theta, par.pos, ref_taper, locs,
                                             x_covariates, smooth.limits,
                                             cholS, z, n, lambda) {
  
  theta_list <- getModelLists(theta = theta, par.pos = par.pos, type = 'diff')
  
  theta_list$std.dev[1] <- 0
  
  ref_taper@entries <- ref_taper@entries * cov_rns_taper_optimized_range(theta = theta_list[-1], 
                                                                         locs = locs, 
                                                                         x_covariates =  x_covariates, 
                                                                         colindices = ref_taper@colindices, 
                                                                         rowpointers = ref_taper@rowpointers,
                                                                         smooth_limits =  smooth.limits)
  
  resid <- z - c(x_covariates %*% theta_list$mean)
  Sigma_cpp <- spam::update.spam.chol.NgPeyton(cholS, ref_taper)
  logdet <- c(spam::determinant.spam.chol.NgPeyton(Sigma_cpp)$modulus)
  
  # check
  return(n * log(2 * pi / n) + n + 2 * logdet + 
           n * log(sum(resid * spam::solve.spam(Sigma_cpp, resid)))) + 
    n * 2 * lambda * exp(theta_list$scale[1]) * 
    sqrt(((smooth.limits[2]-smooth.limits[1])/ 
            (1 + exp(-theta_list$smooth[1])) + 
            smooth.limits[1]))
  
}

#' GetNeg2loglikelihoodProfile
#' @description compute the negative 2 log likelihood based on theta
#'
#' @usage function(theta)
#' @param theta a vector with parameters values
#' @param par.pos a vector with parameters values
#' @param locs a vector with parameters values
#' @param x_covariates a vector with parameters values
#' @param smooth.limits a vector with parameters values
#' @param n a vector with parameters values
#' @param z a vector with parameters values
#' @param lambda amount of regularization
#' @returns neg2loglikelihood 
#' @examples 
#' rnorm(1)
#' @author Federico Blasi
GetNeg2loglikelihoodProfile <- function(theta, par.pos, locs, x_covariates, 
                                        smooth.limits, n, z, x_betas, lambda) {
  
  theta_list <- svcov::getModelLists(theta = theta, par.pos = par.pos, type = 'diff')
  
  Sigma_cpp <- svcov::cov_rns(theta = theta_list[-1], locs = locs, 
                              x_covariates =  x_covariates,
                              smooth_limits = smooth.limits)
  
  # Compute P
  SigmaX <- solve(Sigma_cpp, x_betas)
  P_tete <- solve(Sigma_cpp) - SigmaX %*% solve(t(x_betas) %*% SigmaX) %*% t(SigmaX)
  
  # Quadratic Term
  quad_sum <- t(z) %*% P_tete %*% z
  
  cholS <- base::chol(Sigma_cpp)
  logdet <- sum(log(diag(cholS)))
  
  temp_thetas <- svcov::getModelLists(theta = theta, par.pos = par.pos, type = 'classic')
  
  return(n * log(2 * pi) + 2 * logdet + quad_sum + 
           2 * n * lambda * exp(theta_list$scale[1]) * 
           sqrt(((smooth.limits[2]-smooth.limits[1])/ 
                       (1 + exp(-theta_list$smooth[1])) + 
                       smooth.limits[1])))
}

#' GetNeg2loglikelihood
#' @description compute the negative 2 log likelihood based on theta
#'
#' @usage function(theta)
#' @param theta a vector with parameters values
#' @param par.pos a vector with parameters values
#' @param locs a vector with parameters values
#' @param x_covariates a vector with parameters values
#' @param smooth.limits a vector with parameters values
#' @param n a vector with parameters values
#' @param z a vector with parameters values
#' @param lambda amount of regularization
#' @returns neg2loglikelihood 
#' @examples 
#' rnorm(1)
#' @author Federico Blasi
GetNeg2loglikelihood <- function(theta, 
                                 par.pos, 
                                 locs, 
                                 x_covariates, 
                                 smooth.limits, 
                                 n, 
                                 z,
                                 lambda) {
  
  sum_logliks <- 0
  
  theta_list <- svcov::getModelLists(theta = theta, par.pos = par.pos, type = 'diff')
  
  Sigma_cpp <- svcov::cov_rns(theta = theta_list[-1], locs = locs, x_covariates =  x_covariates,
                              smooth_limits = smooth.limits)
  
  possibleError <- tryCatch(cholS <- base::chol(Sigma_cpp), error = function(e) e)
  if (inherits(possibleError, "error")) {
    return(1e+06)
  } else{
    for(ii in 1:dim(z)[2]){
      
      resid <- c(z[,ii]) - c(x_covariates %*% theta_list$mean)
      
      logdet <- sum(log(diag(cholS)))
      
      sum_logliks <- sum_logliks + n * log(2 * pi) + 2 * logdet + sum(resid * backsolve(cholS,
                                                                                        forwardsolve(cholS, 
                                                                                                     resid, 
                                                                                                     transpose = TRUE, 
                                                                                                     upper.tri = TRUE),
                                                                                        n)) 
      
    }
    
    return(sum_logliks + 2 * dim(z)[2] * n * lambda * exp(theta_list$scale[1]) * 
             sqrt(((smooth.limits[2] - smooth.limits[1])/ 
                         (1 + exp(-theta_list$smooth[1])) + 
                         smooth.limits[1]))
           )
  }
}
