
#' Optimizer for coco objects
#' 
#' @description 
#' Estimation the spatial model parameters using the L-BFGS-B optimizer \[1\].
#' 
#' @usage 
#' cocoOptim(coco.object, boundaries = list(), ncores = "auto", safe = TRUE,
#' optim.type, optim.control)
#' 
#' @param coco.object (\code{S4}) A \link{coco} object.
#' @param boundaries (\code{list}) If provided, a list containing lower, initial, and upper values for the parameters, as defined by \link{getBoundaries}. If not provided, these values are automatically computed with global lower and upper bounds set to -2 and 2.
#' @param ncores (\code{character} or \code{integer}) The number of threads to use for the optimization. If set to `"auto"`, the number of threads is chosen based on system capabilities or a fraction of the available cores.
#' @param optim.type (\code{character}) The optimization approach. Options include:
#' @param safe (\code{logical}) If `TRUE`, the function avoids Cholesky decomposition errors due to ill-posed covariance matrices by returning a pre-defined large value. Defaults to `TRUE`.
#' \itemize{
#'   \item \code{"ml"}: Classical Maximum Likelihood estimation.
#'   \item \code{"pml"}: Profile Maximum Likelihood, factoring out the spatial trend for dense objects or the global marginal variance parameter for sparse objects.
#' }
#' @param optim.control (\code{list}) A list of settings to be passed to the \link[optimParallel]{optimParallel} function \[2\].
#' @returns (\code{S4}) An optimized S4 object of class \code{coco}.
#' @author Federico Blasi
#' @seealso [\link[optimParallel]{optimParallel}]
#' @references 
#' \[1\] Byrd, Richard H., et al. \emph{"A limited memory algorithm for bound constrained optimization."} 
#' SIAM Journal on scientific computing 16.5 (1995): 1190-1208.
#' 
#' \[2\] Gerber, Florian, and Reinhard Furrer. \emph{"optimParallel: An R package providing a parallel version of the L-BFGS-B optimization method."} 
#' R Journal 11.1 (2019): 352-358.
#' @examples
#' \dontrun{
#' model.list <- list('mean' = 0,
#'                    'std.dev' = formula( ~ 1 + cov_x + cov_y),
#'                    'scale' = formula( ~ 1 + cov_x + cov_y),
#'                    'aniso' = 0,
#'                    'tilt' = 0,
#'                    'smooth' = 3/2,
#'                    'nugget' = -Inf)
#'                    
#' coco_object <- coco(type = 'dense',
#'                     data = holes[[1]][1:100,],
#'                     locs = as.matrix(holes[[1]][1:100,1:2]),
#'                     z = holes[[1]][1:100,]$z,
#'                     model.list = model.list)
#'                     
#' optim_coco <- cocoOptim(coco_object,
#' boundaries = getBoundaries(coco_object,
#' lower.value = -3, 3))
#' 
#' plotOptimInfo(optim_coco)
#' 
#' plot(optim_coco)
#' 
#' plot(optim_coco, type = 'ellipse')
#' 
#' plot(optim_coco, type = 'correlations', index = c(2,3,5))
#' 
#' summary(optim_coco)
#'  
#' getEstims(optim_coco)
#' 
#' }
#' 
cocoOptim <- function(coco.object, boundaries = list(), 
                      ncores = "auto", safe = TRUE,
                      optim.type = "ml",
                      optim.control = NULL){
  
  # Init objects
  if(T){
    
    # if optim.control not provided, then some general Optim.control is provided
    optim.control <- if (is.null(optim.control)) {
      getOption("cocons.Optim.Control")
    } else {
      .cocons.update.optim.control(optim.control)
    }
    
    if(is.character(ncores)){
      ncores <- .cocons.set.ncores(coco.object, optim.control)
    } else{
      .cocons.check.ncores(ncores)
    }
    
    cat("using ",ncores," threads.\n")
    
    designMatrix <- cocons::getDesignMatrix(
      model.list = coco.object@model.list,
      data = coco.object@data
    )
    
    # Check boundaries
    if(length(boundaries) == 0){
      boundaries <- cocons::getBoundaries(
        x = coco.object, 
        lower.value = -2,
        upper.value = 2)
    } else{
      .cocons.check.boundaries(coco.object, boundaries)
    }
    
    # Skip scaling
    if(any(colnames(coco.object@data[,coco.object@info$skip.scale]) %in% 
           colnames(coco.object@data))){
      tmp_info <- .cocons.setDesignMatrixCat(coco.object, designMatrix = designMatrix)
      tmp_values <- tmp_info$tmp_values
      mod_DM <- tmp_info$mod_DM
    } else {
      tmp_values <- cocons::getScale(designMatrix$model.matrix)
      mod_DM <- tmp_values$std.covs      
    }
    
    optim.type <- tolower(optim.type)
    
    # Suggested by OptimParallel
    if(tolower(.Platform$OS.type) != "windows"){
      cl <- parallel::makeCluster(spec = ncores, type = "FORK", outfile = "")
    } else{
      cl <- parallel::makeCluster(spec = ncores, outfile = "")
    }
    
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
  }
  
  if (coco.object@type == "dense") {
    
    if(any(c(coco.object@info$lambda.betas, coco.object@info$lambda.Sigma) > 0) && optim.type == "ml"){
      
      if(is.null(coco.object@info$sparse.point)){
        coco.object@info$sparse.point <- getOption("cocons.sparse.point")
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
        "smooth.limits" = coco.object@info$smooth.limits,
        "z" = coco.object@z,
        "x_covariates" = mod_DM,
        "par.pos" = designMatrix$par.pos,
        "lambda" = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg),
        "locs" = coco.object@locs,
        "safe" = safe
      )
      
      rm(designMatrix)
      
      # Call optim routine
      output_dense <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      coco_pen <- cocons:::.cocons.update.coco.first.step(coco.object, output_dense, boundaries)
      
      # .cocons.check.convergence(output_dense, boundaries)
      
      # Create again design matrix
      
      if(T){
        
        designMatrix <- cocons::getDesignMatrix(
          model.list = coco_pen@model.list,
          data = coco_pen@data
        )
        
        # Skip scaling
        if(any(colnames(coco_pen@data[,coco_pen@info$skip.scale]) %in% 
               colnames(coco_pen@data))){
          tmp_info <- .cocons.setDesignMatrixCat(coco_pen, designMatrix = designMatrix)
          tmp_values <- tmp_info$tmp_values
          mod_DM <- tmp_info$mod_DM
        } else {
          tmp_values <- cocons::getScale(designMatrix$model.matrix)
          mod_DM <- tmp_values$std.covs      
        }
        
        args_optim <- list(
          "fn" = cocons::GetNeg2loglikelihood,
          "method" = "L-BFGS-B",
          "lower" = coco_pen@info$boundaries$theta_lower,
          "par" = coco_pen@info$boundaries$theta_init,
          "upper" = coco_pen@info$boundaries$theta_upper,
          "n" = dim(coco_pen@z)[1],
          "smooth.limits" = coco_pen@info$smooth.limits,
          "z" = coco_pen@z,
          "x_covariates" = mod_DM,
          "par.pos" = designMatrix$par.pos,
          "lambda" = c(0, 0, coco_pen@info$lambda.reg),
          "locs" = coco_pen@locs,
          "safe" = safe
        )  
      }
      
      # 
      
      output_dense_two <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      coco_pen@output <- output_dense_two
      
      boundaries_two <- list(theta_init = args_optim$par,
                             theta_upper = args_optim$upper,
                             theta_lower = args_optim$lower)
      
      .cocons.check.convergence(output_dense_two, boundaries_two)
      
      # END ADDED
      
      coco_pen@info$boundaries <- boundaries_two
      coco_pen@info$mean.vector <- tmp_values$mean.vector
      coco_pen@info$sd.vector <- tmp_values$sd.vector
      coco_pen@info$optim.type <- "ml"
      coco_pen@info$safe <- safe
      coco_pen@info$call <- match.call()
      
      return(coco_pen)
      
    }
      
    if(optim.type == "ml"){

      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("cocons"))
      
      args_optim <- list(
        "fn" = cocons::GetNeg2loglikelihood,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "n" = dim(coco.object@z)[1],
        "smooth.limits" = coco.object@info$smooth.limits,
        "z" = coco.object@z,
        "x_covariates" = mod_DM,
        "par.pos" = designMatrix$par.pos,
        "lambda" = c(0, 0, coco.object@info$lambda.reg),
        "locs" = coco.object@locs,
        "safe" = safe
      )
      
      rm(designMatrix)

      # Call optim routine
      output_dense <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      .cocons.check.convergence(output_dense, boundaries)
      
      coco.object@output <- output_dense
      
      coco.object@info$boundaries <- boundaries
      coco.object@info$mean.vector <- tmp_values$mean.vector
      coco.object@info$sd.vector <- tmp_values$sd.vector
      coco.object@info$optim.type <- "ml"
      coco.object@info$safe <- safe
      coco.object@info$call <- match.call()
      
      return(coco.object)
      
    }
    
    if(optim.type == "pml" || optim.type == "reml"){
      
      if(!is.logical(designMatrix$par.pos$mean)){stop("Profile ML or Restricted ML only available when considering covariates in the mean.")}

      x_betas <- mod_DM[, designMatrix$par.pos$mean]
      
      tmp_boundaries <- boundaries
      
      # profiling betas
      boundaries$theta_init <- boundaries$theta_init[-c(1:(sum(designMatrix$par.pos$mean)))]
      boundaries$theta_upper <- boundaries$theta_upper[-c(1:(sum(designMatrix$par.pos$mean)))]
      boundaries$theta_lower <- boundaries$theta_lower[-c(1:(sum(designMatrix$par.pos$mean)))]
      
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("cocons"))
      
      # factoring out betas
      if(T){
        tmp_par_pos <- designMatrix$par.pos
        tmp_par_pos$mean <- rep(FALSE, length(tmp_par_pos$mean))    
      }

      args_optim <- list(
        "fn" = switch(optim.type,
                      pml = cocons::GetNeg2loglikelihoodProfile,
                      reml = cocons::GetNeg2loglikelihoodREML),
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "n" = dim(coco.object@z)[1],
        "smooth.limits" = coco.object@info$smooth.limits,
        "z" = switch(optim.type,
                     pml = coco.object@z,
                     reml = (diag(dim(mod_DM)[1]) - mod_DM %*% solve(crossprod(mod_DM),t(mod_DM))) %*% coco.object@z), 
        "x_covariates" = mod_DM,
        "par.pos" = tmp_par_pos,
        "lambda" = c(0, 0, coco.object@info$lambda.reg),
        "locs" = coco.object@locs,
        "x_betas" = x_betas,
        "safe" = safe
      )
      
      # Call optim routine
      output_dense <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      theta_list <- cocons::getModelLists(theta = output_dense$par, 
                                         par.pos = args_optim$par.pos, type = "diff")
      
      Sigma_cpp <- cocons::cov_rns(theta = theta_list[-1], locs = args_optim$locs,
                                  x_covariates  =  args_optim$x_covariates,
                                  smooth_limits = args_optim$smooth.limits)
      
      # Compute Betas

      if(T){
        L <- chol(Sigma_cpp)
        V <- backsolve(L, forwardsolve(L, x_betas,
                                       transpose = TRUE, 
                                       upper.tri = TRUE))
        W <- crossprod(x_betas, V)
        betass <- c(solve(W, t(V)) %*% rowSums(coco.object@z)) / dim(coco.object@z)[2]
        names(betass) <- colnames(x_betas)
      }
      
      output_dense$par <- c(betass, output_dense$par)
      
      tmp_boundaries$theta_upper[1:(sum(designMatrix$par.pos$mean))] <- rep(Inf,sum(designMatrix$par.pos$mean))
      tmp_boundaries$theta_lower[1:(sum(designMatrix$par.pos$mean))] <- rep(-Inf,sum(designMatrix$par.pos$mean))
      
      .cocons.check.convergence(output_dense, tmp_boundaries)
      
      coco.object@output <- output_dense
      coco.object@info$boundaries <- tmp_boundaries
      coco.object@info$mean.vector <- tmp_values$mean.vector
      coco.object@info$sd.vector <- tmp_values$sd.vector
      coco.object@info$optim.type <- switch(optim.type, pml = "pml", reml = "reml")
      coco.object@info$safe <- safe
      coco.object@info$call <- match.call()
      
      return(coco.object)
      
    }
    
  }
  
  if (coco.object@type == "sparse") {
    
    if(any(c(coco.object@info$lambda.betas, coco.object@info$lambda.Sigma) > 0) && optim.type == "ml"){
      
      if(is.null(coco.object@info$sparse.point)){
        coco.object@info$sparse.point <- getOption("cocons.sparse.point")
      }
      
      # taper
      ref_taper <- coco.object@info$taper(
        spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL),
        theta = c(coco.object@info$delta, 1)
      )
      
      print(summary(ref_taper))
      
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("cocons"))
      parallel::clusterEvalQ(cl, library("spam"))
      #parallel::clusterEvalQ(cl, library("spam64"))
      parallel::clusterEvalQ(cl, options(spam.cholupdatesingular = "error"))
      parallel::clusterEvalQ(cl, options(spam.cholsymmetrycheck = FALSE))
      
      args_optim <- list(
        "fn" = cocons::GetNeg2loglikelihoodTaper,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "par.pos" = designMatrix$par.pos,
        "ref_taper" = ref_taper,
        "locs" = coco.object@locs,
        "x_covariates" = mod_DM,
        "smooth.limits" = coco.object@info$smooth.limits,
        "cholS" = spam::chol.spam(ref_taper),
        "z" = coco.object@z,
        "n" = dim(coco.object@z)[1],
        "lambda" = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg),
        "safe" = safe
      )
      
      # Call optim routine
      output_taper <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      coco_pen <- .cocons.update.coco.first.step(coco.object,output_taper, boundaries)
      
      # Create again design matrix
      
      if(T){
        
        designMatrix <- cocons::getDesignMatrix(
          model.list = coco_pen@model.list,
          data = coco_pen@data
        )
        
        # Skip scaling
        if(any(colnames(coco_pen@data[,coco_pen@info$skip.scale]) %in% 
               colnames(coco_pen@data))){
          tmp_info <- .cocons.setDesignMatrixCat(coco_pen, designMatrix = designMatrix)
          tmp_values <- tmp_info$tmp_values
          mod_DM <- tmp_info$mod_DM
        } else {
          tmp_values <- cocons::getScale(designMatrix$model.matrix)
          mod_DM <- tmp_values$std.covs      
        }
        
        args_optim <- list(
          "fn" = cocons::GetNeg2loglikelihoodTaper,
          "method" = "L-BFGS-B",
          "lower" = coco_pen@info$boundaries$theta_lower,
          "par" = coco_pen@info$boundaries$theta_init,
          "upper" = coco_pen@info$boundaries$theta_upper,
          "par.pos" = designMatrix$par.pos,
          "ref_taper" = ref_taper,
          "locs" = coco_pen@locs,
          "x_covariates" = mod_DM,
          "smooth.limits" = coco_pen@info$smooth.limits,
          "cholS" = spam::chol.spam(ref_taper),
          "z" = coco_pen@z,
          "n" = dim(coco.object@z)[1],
          "lambda" = c(0, 0, coco_pen@info$lambda.reg),
          "safe" = safe
        )
        
      }
        
      output_taper_two <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      coco_pen@output <- output_taper_two
      
      boundaries_two <- list(theta_init = args_optim$par,
                             theta_upper = args_optim$upper,
                             theta_lower = args_optim$lower)
      
      .cocons.check.convergence(output_taper_two, boundaries_two)
      
      coco_pen@info$boundaries <- boundaries_two
      coco_pen@info$mean.vector <- tmp_values$mean.vector
      coco_pen@info$sd.vector <- tmp_values$sd.vector
      coco_pen@info$optim.type <- "ml"
      coco_pen@info$safe <- safe
      coco_pen@info$call <- match.call()
      
      return(coco_pen)
        
      }
      
    if(optim.type == "ml"){
      
      # taper
      ref_taper <- coco.object@info$taper(
        spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL),
        theta = c(coco.object@info$delta, 1)
      )
      
      print(summary(ref_taper))

      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("cocons"))
      parallel::clusterEvalQ(cl, library("spam"))
      #parallel::clusterEvalQ(cl, library("spam64"))
      parallel::clusterEvalQ(cl, options(spam.cholupdatesingular = "error"))
      parallel::clusterEvalQ(cl, options(spam.cholsymmetrycheck = FALSE))
      
      args_optim <- list(
        "fn" = cocons::GetNeg2loglikelihoodTaper,
        "method" = "L-BFGS-B",
        "lower" = boundaries$theta_lower,
        "par" = boundaries$theta_init,
        "upper" = boundaries$theta_upper,
        "par.pos" = designMatrix$par.pos,
        "ref_taper" = ref_taper,
        "locs" = coco.object@locs,
        "x_covariates" = mod_DM,
        "smooth.limits" = coco.object@info$smooth.limits,
        "cholS" = spam::chol.spam(ref_taper),
        "z" = coco.object@z,
        "n" = dim(coco.object@z)[1],
        "lambda" = c(0, 0, coco.object@info$lambda.reg),
        "safe" = safe
      )
      
      # Call optim routine
      output_taper <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      .cocons.check.convergence(output_taper, boundaries)
      
      coco.object@output <- output_taper
      coco.object@info$boundaries <- boundaries
      coco.object@info$mean.vector <- tmp_values$mean.vector
      coco.object@info$sd.vector <- tmp_values$sd.vector
      coco.object@info$optim.type <- "ml"
      coco.object@info$safe <- safe
      coco.object@info$call <- match.call()
      
      return(coco.object)
    }
    
    if(optim.type == "pml"){
      
      if(!is.logical(designMatrix$par.pos$std.dev)){stop("at least a global sigma needs to be estimated for sparse pml coco objects.")}
      
      # taper
      ref_taper <- coco.object@info$taper(
        spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL),
        theta = c(coco.object@info$delta, 1)
      )
      
      print(summary(ref_taper))

      parallel::setDefaultCluster(cl = cl)
      parallel::clusterEvalQ(cl, library("cocons"))
      parallel::clusterEvalQ(cl, library("spam"))
      #parallel::clusterEvalQ(cl, library("spam64"))
      parallel::clusterEvalQ(cl, options(spam.cholupdatesingular = "error"))
      parallel::clusterEvalQ(cl, options(spam.cholsymmetrycheck = FALSE))
      
      tmp_par.pos <- designMatrix$par.pos
      
      tmp_par.pos$std.dev[1] <- FALSE

      boundaries_temp <- boundaries
      
      first_sigma <- if(is.logical(designMatrix$par.pos$mean)){
        first_sigma <- sum(designMatrix$par.pos$mean) + 1
      } else{ first_sigma <- 1}
      
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
        "x_covariates" = mod_DM,
        "smooth.limits" = coco.object@info$smooth.limits,
        "cholS" = spam::chol.spam(ref_taper),
        "z" = coco.object@z,
        "n" = dim(coco.object@z)[1],
        "lambda" = c(0, 0, coco.object@info$lambda.reg),
        "safe" = safe
      )
      
      # Call optim routine
      output_taper <- do.call(what = optimParallel::optimParallel, args = c(
        args_optim,
        optim.control
      ))
      
      if(T){
        
        theta_list <- getModelLists(theta = output_taper$par, par.pos = args_optim$par.pos, type = "diff")
        
        args_optim$ref_taper@entries <- args_optim$ref_taper@entries * cov_rns_taper(theta = theta_list[-1], 
                                                                                                     locs = args_optim$locs, 
                                                                                                     x_covariates =  args_optim$x_covariates, 
                                                                                                     colindices = args_optim$ref_taper@colindices, 
                                                                                                     rowpointers = args_optim$ref_taper@rowpointers,
                                                                                                     smooth_limits =  args_optim$smooth.limits)
        
        cum_qprod <- sum(apply(args_optim$z, MARGIN = 2, FUN = function(x){
          resid <- c(x - args_optim$x_covariates %*% theta_list$mean)
          SigmaX <- spam::solve(args_optim$ref_taper, resid)
          c(resid %*% SigmaX)
        }))
        
        sigma_0 <- cum_qprod / (args_optim$n * dim(args_optim$z)[2])
        
      }
      
      first_scale <- if(is.logical(designMatrix$par.pos$mean)){
        sum(args_optim$par.pos$mean) + sum(args_optim$par.pos$std.dev) + 1
      } else{
        sum(args_optim$par.pos$std.dev) + 1
      }
      
      # what if scale is fixed (?) -- implemented v0.1
        
      tmp_global_scale <- if(is.logical(args_optim$par.pos$scale)){output_taper$par[first_scale]} else{
        tmp_global_scale <- args_optim$par.pos$scale[1]
      }
      
      first_par <- log(sigma_0) + tmp_global_scale
      names(first_par) <- "std.dev.limits"
      second_par <- log(sigma_0) - tmp_global_scale
      names(second_par) <- "scale.limits"
        
      if(is.logical(args_optim$par.pos$mean)){
        new_pars <- output_taper$par[1:(sum(args_optim$par.pos$mean))]
      } else{new_pars <- numeric(0)}
      
      new_pars <- c(new_pars, first_par)
      
      if(is.logical(args_optim$par.pos$std.dev)){
        new_pars <- c(new_pars, output_taper$par[first_sigma:(first_sigma + sum(args_optim$par.pos$std.dev) - 1)])
      }
      
      if(is.logical(args_optim$par.pos$scale)){
        new_pars <- c(new_pars, second_par, output_taper$par[(first_scale+1):(first_scale + sum(args_optim$par.pos$scale) - 1)])
        first_smooth <- first_scale + if(is.logical(args_optim$par.pos$scale)){sum(args_optim$par.pos$scale)}
      } else{
        first_smooth <- first_scale
      }
      
      if(is.logical(args_optim$par.pos$smooth)){
        new_pars <- c(new_pars, output_taper$par[first_smooth:(first_smooth + sum(args_optim$par.pos$smooth) - 1)])
        first_nugget <- first_smooth + if(is.logical(args_optim$par.pos$smooth)){sum(args_optim$par.pos$smooth)}
      } else{first_nugget <- first_smooth}
      
      if(is.logical(args_optim$par.pos$nugget)){
        new_pars <- c(new_pars, output_taper$par[first_nugget:(first_nugget + sum(args_optim$par.pos$nugget) - 1)])
      }
      
      output_taper$par <- new_pars
      
      # fix names output_taper
      
      boundaries_temp$theta_lower[first_sigma] <- -Inf
      boundaries_temp$theta_upper[first_sigma] <- Inf
      
      if(is.logical(args_optim$par.pos$scale)){
        tmp_value <- boundaries_temp$theta_upper[first_scale]
        boundaries_temp$theta_upper[first_scale + 1] <- log(sigma_0) - boundaries_temp$theta_lower[first_scale]
        boundaries_temp$theta_lower[first_scale + 1] <- log(sigma_0) - tmp_value        
      }
      
      .cocons.check.convergence(output_taper, boundaries_temp)
      
      coco.object@output <- output_taper
      coco.object@info$boundaries <- boundaries_temp
      coco.object@info$mean.vector <- tmp_values$mean.vector
      coco.object@info$sd.vector <- tmp_values$sd.vector
      coco.object@info$optim.type <- "pml"
      coco.object@info$safe <- safe
      coco.object@info$call <- match.call()
      
      return(coco.object)
    }
    
  }
  
}
  