
#' Optimizer for object class svcov
#' @description Optimizer based on Optimparallel L-BFGS-B optimizier for svcov class. 
#' @usage svcovOptim(svcov.object, boundaries = list(), ncores = parallel::detectCores())
#' @param svcov.object a svcov object. See ?svcov()
#' @param boundaries if provided, a list with lower, init, and upper values. 
#' if not is computed based on generic fast_init_boundaries()
#' @param ncores number of cores for the optimization routine.
#' @param optim.control number of cores for the optimization routine.
#' @param optim.type Optimization approach.
#' @returns a svcov object with an updated output slot, with extra information 
#' with boundaries information
#' @author Federico Blasi
#' 
svcovOptim <- function(svcov.object, boundaries = list(), 
                       ncores = parallel::detectCores(), 
                       optim.control = NULL, optim.type = 'mle', ...){
  
  if(length(boundaries) > 0){
    .svcov.check.boundaries(boundaries)
  }
  
  .svcov.check.ncores(ncores)
  .svcov.check.z(svcov.object@z) # again checking here due to object might be used to simulate and then attached the realization to x@z
  
  if (svcov.object@type == "dense") {
    
    if(optim.type == 'mle'){
      
      designMatrix <- svcov::getDesignMatrix(
        model.list = svcov.object@model.list,
        data = svcov.object@data
      )
      
      # If boundaries not provided, then some general boundaries are set
      if (length(boundaries) == 0) {
        boundaries <- svcov::getBoundaries(
          x = svcov.object, 
          lower.value = -2,
          upper.value = 2
        )
      }
      
      if(!is.null(svcov.object@info$cat.vars)){
        
        to_not_std <- colnames(svcov.object@data)[svcov.object@info$cat.vars]
        to_avoid_std <- colnames(designMatrix$model.matrix) %in% to_not_std
        
        tmp_values <- svcov::getScale(designMatrix$model.matrix[,!to_avoid_std])
        
        empty_matrix <- matrix(0, ncol = dim(designMatrix$model.matrix)[2], 
                               nrow = dim(designMatrix$model.matrix)[1])
        
        empty_matrix[, !to_avoid_std] <- tmp_values$std.covs
        empty_matrix[, to_avoid_std] <- designMatrix$model.matrix[, to_avoid_std]
        
        mean_vector_empty <- rep(0,dim(designMatrix$model.matrix)[2])
        mean_sd_empty <- rep(1,dim(designMatrix$model.matrix)[2])
        
        mean_vector_empty[!to_avoid_std] <- tmp_values$mean.vector
        mean_sd_empty[!to_avoid_std] <- tmp_values$sd.vector
        
        tmp_values$mean.vector <- mean_vector_empty
        tmp_values$sd.vector <- mean_sd_empty
        
      } else{
        tmp_values <- svcov::getScale(designMatrix$model.matrix)
        empty_matrix <- tmp_values$std.covs
      }
      
      # Suggested by OptimParallel
      if(tolower(.Platform$OS.type) != "windows"){
        cl <- parallel::makeCluster(spec = ncores, type = "FORK", outfile = "")
      } else{
        cl <- parallel::makeCluster(spec = ncores, outfile = "")
      }
      
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("svcov"))
      
      args_optim <- list(
        "fn" = svcov::GetNeg2loglikelihood,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "n" = dim(svcov.object@z)[1],
        "smooth.limits" = svcov.object@info$smooth_limits,
        "z" = svcov.object@z,
        "x_covariates" = empty_matrix,
        "par.pos" = designMatrix$par.pos,
        "lambda" = svcov.object@info$lambda,
        "locs" = svcov.object@locs
      )
      
      rm(designMatrix)
      
      # if optim.control not provided, then some general Optim.control is provided
      if (is.null(optim.control)) {
        optim.control <- getOption("svcov.Optim.Control")
      } else{
        optim.control <- .svcov.update.optim.control(optim.control)
      }
      
      output_dense <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      parallel::stopCluster(cl)
      
      if (any(boundaries$theta_upper == output_dense$par) ||
          any(boundaries$theta_lower == output_dense$par)) {
        warning("at least one of the estimates at the 
              boundaries.")
      }
      
      # Add some warning related to the convergence of the optim routine
      
      svcov.object@output <- output_dense
      svcov.object@info$boundaries <- boundaries
      svcov.object@info$mean.vector <- tmp_values$mean.vector
      svcov.object@info$sd.vector <- tmp_values$sd.vector
      svcov.object@info$optim.type <- 'mle'
      
      return(svcov.object)
      
    }
    
    if(optim.type == 'pmle'){
      
      designMatrix <- svcov::getDesignMatrix(
        model.list = svcov.object@model.list,
        data = svcov.object@data
      )
      
      if(!is.logical(designMatrix$par.pos$mean)){stop('profile ML only available when considering covariates in the trend')}
      
      if (length(boundaries) == 0) {
        boundaries <- svcov::getBoundaries(
          x = svcov.object, 
          lower.value = -3,
          upper.value = 3
        )
      }
      
      # Categorical variables
      # problem here: i can specify categorical values in data but then
      # not use any in formula
      if(!is.null(svcov.object@info$cat.vars)){
        
        to_not_std <- colnames(svcov.object@data)[svcov.object@info$cat.vars]
        to_avoid_std <- colnames(designMatrix$model.matrix) %in% to_not_std
        
        tmp_values <- svcov::getScale(designMatrix$model.matrix[,!to_avoid_std])
        
        empty_matrix <- matrix(0, ncol = dim(designMatrix$model.matrix)[2], 
                               nrow = dim(designMatrix$model.matrix)[1])
        
        empty_matrix[, !to_avoid_std] <- tmp_values$std.covs
        empty_matrix[, to_avoid_std] <- designMatrix$model.matrix[, to_avoid_std]
        
        mean_vector_empty <- rep(0,dim(designMatrix$model.matrix)[2])
        mean_sd_empty <- rep(1,dim(designMatrix$model.matrix)[2])
        
        mean_vector_empty[!to_avoid_std] <- tmp_values$mean.vector
        mean_sd_empty[!to_avoid_std] <- tmp_values$sd.vector
        
        tmp_values$mean.vector <- mean_vector_empty
        tmp_values$sd.vector <- mean_sd_empty
        
      } else{
        tmp_values <- svcov::getScale(designMatrix$model.matrix)
        empty_matrix <- tmp_values$std.covs
      }
      
      if(is.logical(designMatrix$par.pos$mean)){
        x_betas <- getScale(svcov.object)$std.covs[, designMatrix$par.pos$mean] # does not match when categorical variables
      }
      
      tmp_boundaries <- boundaries
      
      # profiling both betas
      boundaries$theta_init <- boundaries$theta_init[-c(1:(sum(designMatrix$par.pos$mean)))]
      boundaries$theta_upper <- boundaries$theta_upper[-c(1:(sum(designMatrix$par.pos$mean)))]
      boundaries$theta_lower <- boundaries$theta_lower[-c(1:(sum(designMatrix$par.pos$mean)))]
      
      # Suggested by OptimParallel
      if(tolower(.Platform$OS.type) != "windows"){
        cl <- parallel::makeCluster(spec = ncores, type = "FORK", outfile = "")
      } else{
        cl <- parallel::makeCluster(spec = ncores, outfile = "")
      }
      
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("svcov"))
      
      # factoring out betas
      tmp_par_pos <- designMatrix$par.pos
      
      if(is.logical(tmp_par_pos$mean)){
        # factoring out all betas
        tmp_par_pos$mean <- rep(FALSE, length(tmp_par_pos$mean))
      }
      
      args_optim <- list(
        "fn" = svcov::GetNeg2loglikelihoodProfile,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "n" = length(svcov.object@z),
        "smooth.limits" = svcov.object@info$smooth_limits,
        "z" = svcov.object@z,
        "x_covariates" = empty_matrix,
        "par.pos" = tmp_par_pos,
        "lambda" = svcov.object@info$lambda,
        "locs" = svcov.object@locs,
        "x_betas" = x_betas
      )
      
      # rm(designMatrix)
      
      # if optim.control not provided, then some general Optim.control is provided
      if (is.null(optim.control)) {
        optim.control <- getOption("svcov.Optim.Control")
      } else{
        optim.control <- .svcov.update.optim.control(optim.control)
      }
      
      output_dense <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      parallel::stopCluster(cl)
      
      theta_list <- svcov::getModelLists(theta = output_dense$par, 
                                         par.pos = args_optim$par.pos, type = 'diff')
      
      Sigma_cpp <- svcov::cov_rns(theta = theta_list[-1], locs = args_optim$locs,
                                  x_covariates  =  args_optim$x_covariates,
                                  smooth_limits = args_optim$smooth.limits)
      
      Sigma_X <- solve(Sigma_cpp, args_optim$x_betas)
      betass <- solve(t(args_optim$x_betas) %*% Sigma_X) %*% t(Sigma_X) %*% args_optim$z
      output_dense$par <- c(betass, output_dense$par)
      
      if (any(tmp_boundaries$theta_upper == output_dense$par, na.rm = T) ||
          any(tmp_boundaries$theta_lower == output_dense$par, na.rm = T)) {
        warning("at least one of the estimates at the 
              boundaries.")
      }
      
      # Add some warning related to the convergence of the optim routine
      
      svcov.object@output <- output_dense
      svcov.object@info$boundaries <- tmp_boundaries
      svcov.object@info$mean.vector <- tmp_values$mean.vector
      svcov.object@info$sd.vector <- tmp_values$sd.vector
      svcov.object@info$optim.type <- 'pmle'
      
      return(svcov.object)
      
    }
    
  }
  
  if (svcov.object@type == "sparse") {
    
    if(optim.type == 'mle'){
      
      designMatrix <- svcov::getDesignMatrix(
        model.list = svcov.object@model.list,
        data = svcov.object@data
      )
      
      if (length(boundaries) == 0) {
        boundaries <- svcov::getBoundaries(
          x = svcov.object, 
          lower.value = -3,
          upper.value = 3
        )
      }
      
      # Categorical variables
      # problem here: i can specify categorical values in data but then
      # not use any in formula
      if(!is.null(svcov.object@info$cat.vars)){
        
        to_not_std <- colnames(svcov.object@data)[svcov.object@info$cat.vars]
        to_avoid_std <- colnames(designMatrix$model.matrix) %in% to_not_std
        
        tmp_values <- svcov::getScale(designMatrix$model.matrix[,!to_avoid_std])
        
        empty_matrix <- matrix(0, ncol = dim(designMatrix$model.matrix)[2], 
                               nrow = dim(designMatrix$model.matrix)[1])
        
        empty_matrix[, !to_avoid_std] <- tmp_values$std.covs
        empty_matrix[, to_avoid_std] <- designMatrix$model.matrix[, to_avoid_std]
        
        mean_vector_empty <- rep(0,dim(designMatrix$model.matrix)[2])
        mean_sd_empty <- rep(1,dim(designMatrix$model.matrix)[2])
        
        mean_vector_empty[!to_avoid_std] <- tmp_values$mean.vector
        mean_sd_empty[!to_avoid_std] <- tmp_values$sd.vector
        
        tmp_values$mean.vector <- mean_vector_empty
        tmp_values$sd.vector <- mean_sd_empty
        
      } else{
        tmp_values <- svcov::getScale(designMatrix$model.matrix)
        empty_matrix <- tmp_values$std.covs
      }
      
      # taper
      ref_taper <- svcov.object@info$taper(
        spam::nearest.dist(svcov.object@locs, delta = svcov.object@info$delta, upper = NULL),
        theta = c(svcov.object@info$delta, 1)
      )
      
      print(summary(ref_taper))
      
      if(tolower(.Platform$OS.type) != "windows"){
        cl <- parallel::makeCluster(spec = ncores, type = "FORK", outfile = "",
                                    useXDR = FALSE)
      } else{
        cl <- parallel::makeCluster(spec = ncores, outfile = "", 
                                    useXDR = FALSE)
      }
      
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("svcov"))
      parallel::clusterEvalQ(cl, library("spam"))
      parallel::clusterEvalQ(cl, library("spam64"))
      parallel::clusterEvalQ(cl, options(spam.cholupdatesingular = "error"))
      parallel::clusterEvalQ(cl, options(spam.cholsymmetrycheck = FALSE)) # check whether i want to add this
      
      # options(spam.cholsymmetrycheck = FALSE)
      
      args_optim <- list(
        "fn" = svcov::GetNeg2loglikelihoodTaper,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "par.pos" = designMatrix$par.pos,
        "ref_taper" = ref_taper,
        "locs" = svcov.object@locs,
        "x_covariates" = empty_matrix,
        "smooth.limits" = svcov.object@info$smooth_limits,
        "cholS" = spam::chol.spam(ref_taper),
        "z" = svcov.object@z,
        "n" = length(svcov.object@z),
        "lambda" = svcov.object@info$lambda
      )
      
      #rm(designMatrix)
      #rm(ref_taper)
      
      # if optim.control not provided, then some general Optim.control is provided
      if (is.null(optim.control)) {
        optim.control <- getOption("svcov.Optim.Control")
      } else{
        optim.control <- .svcov.update.optim.control(optim.control)
      }
      
      output_taper <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      parallel::stopCluster(cl)
      
      if (any(boundaries$theta_upper == output_taper$par) ||
          any(boundaries$theta_lower == output_taper$par)) {
        warning("at least one of the estimates at the 
              boundaries.")
      }
      
      svcov.object@output <- output_taper
      svcov.object@info$boundaries <- boundaries
      svcov.object@info$mean.vector <- tmp_values$mean.vector
      svcov.object@info$sd.vector <- tmp_values$sd.vector
      svcov.object@info$optim.type <- 'mle'
      # Add some warning related to the convergence of the optim routine
      
      return(svcov.object)
    }
    
    if(optim.type == 'pmle'){
      
      designMatrix <- svcov::getDesignMatrix(
        model.list = svcov.object@model.list,
        data = svcov.object@data
      )
      
      if (length(boundaries) == 0) {
        boundaries <- svcov::getBoundaries(
          x = svcov.object, 
          lower.value = -3,
          upper.value = 3
        )
      }
      
      # Categorical variables
      # problem here: i can specify categorical values in data but then
      # not use any in formula
      if(!is.null(svcov.object@info$cat.vars)){
        
        to_not_std <- colnames(svcov.object@data)[svcov.object@info$cat.vars]
        to_avoid_std <- colnames(designMatrix$model.matrix) %in% to_not_std
        
        tmp_values <- svcov::getScale(designMatrix$model.matrix[,!to_avoid_std])
        
        empty_matrix <- matrix(0, ncol = dim(designMatrix$model.matrix)[2], 
                               nrow = dim(designMatrix$model.matrix)[1])
        
        empty_matrix[, !to_avoid_std] <- tmp_values$std.covs
        empty_matrix[, to_avoid_std] <- designMatrix$model.matrix[, to_avoid_std]
        
        mean_vector_empty <- rep(0,dim(designMatrix$model.matrix)[2])
        mean_sd_empty <- rep(1,dim(designMatrix$model.matrix)[2])
        
        mean_vector_empty[!to_avoid_std] <- tmp_values$mean.vector
        mean_sd_empty[!to_avoid_std] <- tmp_values$sd.vector
        
        tmp_values$mean.vector <- mean_vector_empty
        tmp_values$sd.vector <- mean_sd_empty
        
      } else{
        tmp_values <- svcov::getScale(designMatrix$model.matrix)
        empty_matrix <- tmp_values$std.covs
      }
      
      # taper
      ref_taper <- svcov.object@info$taper(
        spam::nearest.dist(svcov.object@locs, delta = svcov.object@info$delta, upper = NULL),
        theta = c(svcov.object@info$delta, 1)
      )
      
      print(summary(ref_taper))
      
      # Suggested by OptimParallel
      if(tolower(.Platform$OS.type) != "windows"){
        cl <- parallel::makeCluster(spec = ncores, type = "FORK", outfile = "",
                                    useXDR = FALSE)
      } else{
        cl <- parallel::makeCluster(spec = ncores, outfile = "", 
                                    useXDR = FALSE)
      }
      
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("svcov"))
      parallel::clusterEvalQ(cl, library("spam"))
      parallel::clusterEvalQ(cl, library("spam64"))
      parallel::clusterEvalQ(cl, options(spam.cholupdatesingular = "error"))
      parallel::clusterEvalQ(cl, options(spam.cholsymmetrycheck = FALSE)) # check whether i want to add this
      
      # options(spam.cholsymmetrycheck = FALSE)
      
      tmp_par.pos <- designMatrix$par.pos
      
      if(is.logical(tmp_par.pos$std.dev)){
        tmp_par.pos$std.dev[1] <- FALSE
      }
      
      boundaries_temp <- boundaries
      
      first_sigma <- sum(designMatrix$par.pos$mean) + 1
      
      boundaries$theta_init <- boundaries$theta_init[-first_sigma]
      boundaries$theta_upper <- boundaries$theta_upper[-first_sigma]
      boundaries$theta_lower <- boundaries$theta_lower[-first_sigma]
      
      args_optim <- list(
        "fn" = svcov::GetNeg2loglikelihoodTaperProfile,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "par.pos" = tmp_par.pos,
        "ref_taper" = ref_taper,
        "locs" = svcov.object@locs,
        "x_covariates" = empty_matrix,
        "smooth.limits" = svcov.object@info$smooth_limits,
        "cholS" = spam::chol.spam(ref_taper),
        "z" = svcov.object@z,
        "n" = length(svcov.object@z),
        "lambda" = svcov.object@info$lambda
      )
      
      #rm(designMatrix)
      #rm(ref_taper)
      
      # if optim.control not provided, then some general Optim.control is provided
      if (is.null(optim.control)) {
        optim.control <- getOption("svcov.Optim.Control")
      } else{
        optim.control <- .svcov.update.optim.control(optim.control)
      }
      
      output_taper <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      parallel::stopCluster(cl)
      
      if(T){
        
        theta_list <- getModelLists(theta = output_taper$par, par.pos = args_optim$par.pos, type = 'diff')
        
        args_optim$ref_taper@entries <- args_optim$ref_taper@entries * cov_rns_taper_optimized_range(theta = theta_list[-1], 
                                                                                                     locs = args_optim$locs, 
                                                                                                     x_covariates =  args_optim$x_covariates, 
                                                                                                     colindices = args_optim$ref_taper@colindices, 
                                                                                                     rowpointers = args_optim$ref_taper@rowpointers,
                                                                                                     smooth_limits =  args_optim$smooth.limits)
        
        resid <- c(args_optim$z - args_optim$x_covariates %*% theta_list$mean)
        
        SigmaX <- spam::solve(args_optim$ref_taper, resid)
        sigma_0 <- c(resid %*% SigmaX) / args_optim$n
        
      }
      
      first_scale <- sum(args_optim$par.pos$mean) + sum(args_optim$par.pos$std.dev) + 1
      
      tmp_global_scale <- output_taper$par[first_scale]
      
      first_par <- log(sigma_0) + tmp_global_scale
      names(first_par) <- 'std.dev.limits'
      second_par <- log(sigma_0) - tmp_global_scale
      
      if(is.logical(args_optim$par.pos$mean)){
        new_pars <- output_taper$par[1:(sum(args_optim$par.pos$mean))]
      } else{new_pars <- numeric(0)}
      
      new_pars <- c(new_pars, first_par)
      
      if(is.logical(args_optim$par.pos$std.dev)){
        new_pars <- c(new_pars, output_taper$par[first_sigma:(first_sigma + sum(args_optim$par.pos$std.dev) - 1)])
      }
      
      if(is.logical(args_optim$par.pos$scale)){
        new_pars <- c(new_pars, second_par, output_taper$par[(first_scale+1):(first_scale + sum(args_optim$par.pos$scale) - 1)])
      }
      
      first_smooth <- first_scale + sum(args_optim$par.pos$scale)
      
      if(is.logical(args_optim$par.pos$smooth)){
        new_pars <- c(new_pars, output_taper$par[first_smooth:(first_smooth + sum(args_optim$par.pos$smooth) - 1)])
      }
      
      first_nugget <- first_smooth + sum(args_optim$par.pos$smooth)
      
      if(sum(args_optim$par.pos$nugget) > 0){
        new_pars <- c(new_pars, output_taper$par[first_nugget:(first_nugget + sum(args_optim$par.pos$nugget) - 1)])
      }
      
      output_taper$par <- new_pars
      
      # fix names output_taper
      
      if (any(boundaries_temp$theta_upper == output_taper$par, na.rm = T) ||
          any(boundaries_temp$theta_lower == output_taper$par, na.rm = T)) {
        warning("at least one of the estimates at the 
              boundaries.")
      }
      
      svcov.object@output <- output_taper
      svcov.object@info$boundaries <- boundaries_temp
      svcov.object@info$mean.vector <- tmp_values$mean.vector
      svcov.object@info$sd.vector <- tmp_values$sd.vector
      svcov.object@info$optim.type <- 'pmle'
      # Add some warning related to the convergence of the optim routine
      
      return(svcov.object)
    }
    
  }
  
}
