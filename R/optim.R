
#' Optimizer for object class coco
#' @description Optimizer based on multi-thread Optimparallel L-BFGS-B optimizers.
#' @details
#'  Current implementations only allow single realizations for 'pmle' optim.type. 
#' @usage cocoOptim(coco.object, boundaries = list(), 
#' ncores = parallel::detectCores(), optim.control, optim.type)
#' @param coco.object a coco object. See ?coco()
#' @param boundaries if provided, a list with lower, init, and upper values. 
#' if not is computed based on generic fast_init_boundaries()
#' @param ncores number of cores to be used for the optimization routine.
#' @param optim.control list with settings to be passed to the optimparallel function
#' @param optim.type Optimization approach: whether 'mle' for classical Maximum Likelihood approach, 
#' or 'pmle' to factor out the spatial trend (when handling 'dense' coco objects), or
#' to factor out the global marginal standard deviation parameter (when considering 'sparse' coco objects)
#' @returns a fitted coco object
#' @author Federico Blasi
#' 
cocoOptim <- function(coco.object, boundaries = list(), 
                       ncores = parallel::detectCores(), 
                       optim.control = NULL, optim.type = "mle"){
  
  if(length(boundaries) > 0){
    .cocons.check.boundaries(boundaries)
  }
  
  .cocons.check.ncores(ncores)
  
  # Init objects
  if(T){
    
    designMatrix <- cocons::getDesignMatrix(
      model.list = coco.object@model.list,
      data = coco.object@data
    )
    
    if (length(boundaries) == 0) {
      boundaries <- cocons::getBoundaries(
        x = coco.object, 
        lower.value = -2,
        upper.value = 2
      )
    }
    
    # Categorical variables
    # problem here: i can specify categorical values in data but then not use any in formula
    if(!is.null(coco.object@info$cat.vars)){
      tmp_info <- .cocons.setDesignMatrixCat(coco.object, designMatrix = designMatrix)
      tmp_values <- tmp_info$tmp_values
      empty_matrix <- tmp_info$empty_matrix
    } else{
      tmp_values <- cocons::getScale(designMatrix$model.matrix)
      empty_matrix <- tmp_values$std.covs
    }
    
  }
  
  if (coco.object@type == "dense") {
    
    if(optim.type == "mle"){
      
      # Suggested by OptimParallel
      if(tolower(.Platform$OS.type) != "windows"){
        cl <- parallel::makeCluster(spec = ncores, type = "FORK", outfile = "")
      } else{
        cl <- parallel::makeCluster(spec = ncores, outfile = "")
      }
      
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("cocons"))
      
      args_optim <- list(
        "fn" = cocons::GetNeg2loglikelihood,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "n" = dim(coco.object@z)[1],
        "smooth.limits" = coco.object@info$smooth_limits,
        "z" = coco.object@z,
        "x_covariates" = empty_matrix,
        "par.pos" = designMatrix$par.pos,
        "lambda" = coco.object@info$lambda,
        "locs" = coco.object@locs
      )
      
      rm(designMatrix)
      
      # if optim.control not provided, then some general Optim.control is provided
      if (is.null(optim.control)) {
        optim.control <- getOption("cocons.Optim.Control")
      } else{
        optim.control <- .cocons.update.optim.control(optim.control)
      }
      
      # Call optim routine
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
      
      coco.object@output <- output_dense
      coco.object@info$boundaries <- boundaries
      coco.object@info$mean.vector <- tmp_values$mean.vector
      coco.object@info$sd.vector <- tmp_values$sd.vector
      coco.object@info$optim.type <- "mle"
      
      return(coco.object)
      
    }
    
    if(optim.type == "pmle"){
      
      if(dim(coco.object@z)[2] > 1){stop('profile ML routines only handle single realizations.')}
      
      if(!is.logical(designMatrix$par.pos$mean)){stop("profile ML only available when considering covariates in the trend")}

      if(is.logical(designMatrix$par.pos$mean)){
        x_betas <- empty_matrix[, designMatrix$par.pos$mean]
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
      parallel::clusterEvalQ(cl, library("cocons"))
      
      # factoring out betas
      tmp_par_pos <- designMatrix$par.pos
      
      if(is.logical(tmp_par_pos$mean)){
        # factoring out all betas
        tmp_par_pos$mean <- rep(FALSE, length(tmp_par_pos$mean))
      }
      
      args_optim <- list(
        "fn" = cocons::GetNeg2loglikelihoodProfile,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "n" = length(coco.object@z),
        "smooth.limits" = coco.object@info$smooth_limits,
        "z" = coco.object@z,
        "x_covariates" = empty_matrix,
        "par.pos" = tmp_par_pos,
        "lambda" = coco.object@info$lambda,
        "locs" = coco.object@locs,
        "x_betas" = x_betas
      )
      
      # rm(designMatrix)
      
      # if optim.control not provided, then some general Optim.control is provided
      if (is.null(optim.control)) {
        optim.control <- getOption("cocons.Optim.Control")
      } else{
        optim.control <- .cocons.update.optim.control(optim.control)
      }
      
      # Call optim routine
      output_dense <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      parallel::stopCluster(cl)
      
      theta_list <- cocons::getModelLists(theta = output_dense$par, 
                                         par.pos = args_optim$par.pos, type = "diff")
      
      Sigma_cpp <- cocons::cov_rns(theta = theta_list[-1], locs = args_optim$locs,
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
      
      coco.object@output <- output_dense
      coco.object@info$boundaries <- tmp_boundaries
      coco.object@info$mean.vector <- tmp_values$mean.vector
      coco.object@info$sd.vector <- tmp_values$sd.vector
      coco.object@info$optim.type <- "pmle"
      
      return(coco.object)
      
    }
    
  }
  
  if (coco.object@type == "sparse") {
    
    if(optim.type == "mle"){
      
      # taper
      ref_taper <- coco.object@info$taper(
        spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL),
        theta = c(coco.object@info$delta, 1)
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
      parallel::clusterEvalQ(cl, library("cocons"))
      parallel::clusterEvalQ(cl, library("spam"))
      parallel::clusterEvalQ(cl, library("spam64"))
      parallel::clusterEvalQ(cl, options(spam.cholupdatesingular = "error"))
      parallel::clusterEvalQ(cl, options(spam.cholsymmetrycheck = FALSE)) # check whether i want to add this
      
      # options(spam.cholsymmetrycheck = FALSE)
      
      args_optim <- list(
        "fn" = cocons::GetNeg2loglikelihoodTaper,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "par.pos" = designMatrix$par.pos,
        "ref_taper" = ref_taper,
        "locs" = coco.object@locs,
        "x_covariates" = empty_matrix,
        "smooth.limits" = coco.object@info$smooth_limits,
        "cholS" = spam::chol.spam(ref_taper),
        "z" = coco.object@z,
        "n" = length(coco.object@z),
        "lambda" = coco.object@info$lambda
      )
      
      #rm(designMatrix)
      #rm(ref_taper)
      
      # if optim.control not provided, then some general Optim.control is provided
      if (is.null(optim.control)) {
        optim.control <- getOption("cocons.Optim.Control")
      } else{
        optim.control <- .cocons.update.optim.control(optim.control)
      }
      
      # Call optim routine
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
      
      coco.object@output <- output_taper
      coco.object@info$boundaries <- boundaries
      coco.object@info$mean.vector <- tmp_values$mean.vector
      coco.object@info$sd.vector <- tmp_values$sd.vector
      coco.object@info$optim.type <- "mle"
      # Add some warning related to the convergence of the optim routine
      
      return(coco.object)
    }
    
    if(optim.type == "pmle"){
      
      if(dim(coco.object@z)[2] > 1){stop('profile ML routines only handle single realizations.')}
      
      # taper
      ref_taper <- coco.object@info$taper(
        spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL),
        theta = c(coco.object@info$delta, 1)
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
      parallel::clusterEvalQ(cl, library("cocons"))
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
        "fn" = cocons::GetNeg2loglikelihoodTaperProfile,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "par.pos" = tmp_par.pos,
        "ref_taper" = ref_taper,
        "locs" = coco.object@locs,
        "x_covariates" = empty_matrix,
        "smooth.limits" = coco.object@info$smooth_limits,
        "cholS" = spam::chol.spam(ref_taper),
        "z" = coco.object@z,
        "n" = length(coco.object@z),
        "lambda" = coco.object@info$lambda
      )
      
      #rm(designMatrix)
      #rm(ref_taper)
      
      # if optim.control not provided, then some general Optim.control is provided
      if (is.null(optim.control)) {
        optim.control <- getOption("cocons.Optim.Control")
      } else{
        optim.control <- .cocons.update.optim.control(optim.control)
      }
      
      # Call optim routine
      output_taper <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      parallel::stopCluster(cl)
      
      if(T){
        
        theta_list <- getModelLists(theta = output_taper$par, par.pos = args_optim$par.pos, type = "diff")
        
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
      names(first_par) <- "std.dev.limits"
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
      
      coco.object@output <- output_taper
      coco.object@info$boundaries <- boundaries_temp
      coco.object@info$mean.vector <- tmp_values$mean.vector
      coco.object@info$sd.vector <- tmp_values$sd.vector
      coco.object@info$optim.type <- "pmle"
      # Add some warning related to the convergence of the optim routine
      
      return(coco.object)
    }
    
  }
  
}
