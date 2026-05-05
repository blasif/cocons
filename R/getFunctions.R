#' Covariance matrix for "coco" class
#' @description Compute the covariance matrix of \code{coco.object}.
#'
#' @usage getCovMatrix(coco.object, type = 'global', index = NULL)
#' @param coco.object \code{(S4)} a fitted [coco()] object.
#' @param type \code{(character)} whether \code{'global'} to retrieve the regular covariance matrix, or \code{'local'} to retrieve global covariance.
#' based on the local aspects of a specific location (not implemented yet).
#' @param index \code{(integer)} index to perform local covariance matrix (not implemented yet).
#' @returns (\code{matrix} or \code{S4}) a n x n covariance matrix (for 'dense' coco objects) or a S4 spam object (for 'sparse' coco objects).
#' @author Federico Blasi
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
#' getCovMatrix(optim_coco)
#' 
#' }
#' 
getCovMatrix <- function(coco.object, type = "global", index = NULL){
  
  x_covs <- cocons::getScale(coco.object)$std.covs
  
  par.pos <- getDesignMatrix(coco.object@model.list, data = coco.object@data)$par.pos
  
  theta_list <- cocons::getModelLists(coco.object@output$par,par.pos = par.pos, 
                                      type = "diff")
  
  if(coco.object@type == "dense"){
    
    if(type == "global"){

      return(cocons::cov_rns(theta = theta_list,
                            locs = coco.object@locs,
                            x_covariates = x_covs,
                            smooth_limits = coco.object@info$smooth.limits))
    }
    
    if(type == "local"){
      
    }
    
  }
  
  if(coco.object@type == "sparse"){
    
    if(type == "global"){

      # taper
      ref_taper <- coco.object@info$taper(
        spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL),
        theta = c(coco.object@info$delta, 1)
      )
      
      ref_taper@entries <- ref_taper@entries * cov_rns_taper(theta = theta_list[-1], 
                                                                                 locs = coco.object@locs, 
                                                                                 x_covariates =  x_covs, 
                                                                                 colindices = ref_taper@colindices, 
                                                                                 rowpointers = ref_taper@rowpointers,
                                                                                 smooth_limits =  coco.object@info$smooth.limits)
      
      return(ref_taper)
      
    }
    
    if(type == "local"){
      
      
    }
    
  }
}

#' Based on a set of predictions computes the Log-Score
#' @description Computes the Log-Score \[1\].
#'
#' @usage getLogScore(z.pred, mean.pred, sd.pred)
#' @param z.pred \code{(numeric vector)}.
#' @param mean.pred \code{(numeric vector)}.
#' @param sd.pred \code{(numeric vector)}.
#' @returns (\code{numeric vector}) retrieves Log-Score.
#' @author Federico Blasi
#' @references
#' \[1\] Gneiting, Tilmann, and Adrian E. Raftery. \emph{"Strictly proper scoring rules, prediction, and estimation."} Journal of the American statistical Association 102.477 (2007): 359-378.
getLogScore <- function(z.pred, mean.pred, sd.pred){
  
  return( (log(2 * pi) + ((z.pred - mean.pred) / sd.pred)^2)/2 + log(sd.pred))
  
}

#' Based on a set of predictions computes the Continuous Ranked Probability Score
#' @description Retrieves the Continuous Ranked Probability Score (CRPS) \[1\].
#'
#' @usage getCRPS(z.pred, mean.pred, sd.pred)
#' @param z.pred \code{(numeric vector)}.
#' @param mean.pred \code{(numeric vector)}.
#' @param sd.pred \code{(numeric vector)}.
#' @returns (\code{numeric vector}) retrieves CRPS.
#' @author Federico Blasi
#' @references 
#' \[1\] Gneiting, Tilmann, and Adrian E. Raftery. \emph{"Strictly proper scoring rules, prediction, and estimation."} Journal of the American statistical Association 102.477 (2007): 359-378.
getCRPS <- function(z.pred, mean.pred, sd.pred){
  
  vector_tmp_z <- (mean.pred - z.pred) / sd.pred
  
  return( sd.pred * (vector_tmp_z * (2 * stats::pnorm(vector_tmp_z) - 1) + 
                       2 * stats::dnorm(vector_tmp_z) - 1 / base::sqrt(pi)))
  
}

#' Based on a specific taper scale (delta), retrieves the density of the covariance matrix.
#' @description Based on a specific taper scale (delta), retrieves the density of the covariance matrix.
#'
#' @usage getDensityFromDelta(coco.object, delta)
#' @param coco.object \code{(S4)} a fitted [coco()] object.
#' @param delta \code{(numeric)} a delta taper scale (delta).
#' @returns (\code{numeric vector}) the associate density of the tapered covariance matrix.
#' @author Federico Blasi
getDensityFromDelta <- function(coco.object, delta){
  
  if(coco.object@type != "sparse"){stop("only for sparse coco objects.")}
  
  if(length(coco.object@output) == 0){
    coco.object@output$par <- rep(0, .cocons.get.npars(coco.object))
  }
  
  coco.object@info$delta <- delta
  
  return(summary(getCovMatrix(coco.object))$density / 100)
  
}

#' Evaluates the spatially-varying functions from a coco object at locs
#' @description Evaluates the spatially-varying functions of the nonstationary spatial structure.
#'
#' @usage getSpatEffects(coco.object)
#' @param coco.object \code{(S4)} a fitted coco S4 object.
#' @returns (\code{list}) a list with the different estimated surfaces.
#' @author Federico Blasi
getSpatEffects <- function(coco.object){
  
  tmp_info <- cocons::getDesignMatrix(model.list = coco.object@model.list, 
                                     data = coco.object@data)
  
  theta_list <- cocons::getModelLists(
    theta = coco.object@output$par, par.pos = tmp_info$par.pos,
    type = "diff"
  )
  
  X_std <- cocons::getScale(tmp_info$model.matrix,
                           mean.vector = coco.object@info$mean.vector,
                           sd.vector = coco.object@info$sd.vector
  )
  
  if(coco.object@type == "dense"){
    
    tp_se <- exp(0.5 * X_std$std.covs %*% theta_list$std.dev)
    tp_ga <- exp(X_std$std.covs %*% theta_list$aniso)
    tp_tl <- pi / (1 + exp(-X_std$std.covs %*% theta_list$tilt))
    tp_smooth <- (coco.object@info$smooth.limits[2] - coco.object@info$smooth.limits[1]) / (1 + exp(-X_std$std.covs %*% theta_list$smooth)) + coco.object@info$smooth.limits[1]
    tp_ng <- exp(X_std$std.covs %*% theta_list$nugget)
    tp_mr_x <- sin(tp_tl) * exp(X_std$std.covs %*% theta_list$scale)
    tp_mr_y <- sin(tp_tl) * exp(X_std$std.covs %*% theta_list$scale) * exp(X_std$std.covs %*% theta_list$aniso)

    
    tp_angle <- atan2(x = 2 * tp_ga^(0.5) * cospi(tp_tl / pi), 
          y = tp_ga - 1 + sqrt((tp_ga + 1)^2 - 4 * tp_ga * sinpi(tp_tl / pi)^2)) # * 180 / pi

    return(list("sd" = tp_se,
                "scale_x" = tp_mr_x,
                "scale_y" = tp_mr_y,
                "aniso" = tp_ga,
                "tilt" = tp_tl,
                "angle" = tp_angle,
                "smooth" = tp_smooth,
                "nugget" = tp_ng))
    
  }
  
  if(coco.object@type == "sparse"){
    
    tp_se <- exp(0.5 * X_std$std.covs %*% theta_list$std.dev)
    tp_smooth <- (coco.object@info$smooth.limits[2] - coco.object@info$smooth.limits[1]) / (1 + exp(-X_std$std.covs %*% theta_list$smooth)) + coco.object@info$smooth.limits[1]
    tp_ng <- exp(X_std$std.covs %*% theta_list$nugget)
    tp_mr <- exp(X_std$std.covs %*% theta_list$scale)

    return(list("sd" = tp_se,
                "scale_x" = tp_mr,
                "smooth" = tp_smooth,
                "nugget" = tp_ng))
    
  }

}

#' Computes the spatial mean of a (fitted) coco object
#' @description Computes the spatial mean of the (fitted) coco object.
#'
#' @usage getSpatMean(coco.object)
#' @param coco.object \code{(S4)} a fitted coco S4 object.
#' @returns (\code{numeric vector}) a vector with the adjusted trend.
#' @author Federico Blasi
getSpatMean <- function(coco.object){
  tmp_scaled <- getScale(coco.object)$std.covs
  return(tmp_scaled %*% getEstims(coco.object)$mean)
}

#' Compute approximate confidence intervals for a coco object
#' @description Compute approximate confidence intervals for a (fitted) coco object.
#'
#' @usage getCIs(coco.object, inv.hess, alpha = 0.95)
#' @param coco.object \code{(S4)} a fitted coco S4 object.
#' @param inv.hess \code{(matrix)} Inverse of the Hessian. \link{getHessian}.
#' @param alpha \code{(numeric)} confidence level.
#' @returns (\code{numeric matrix}) a matrix with approximate confidence intervals for each parameter in the model.
#' @author Federico Blasi
getCIs <- function(coco.object, inv.hess, alpha = 0.95){
  
  if(alpha >= 1 || alpha <= 0){stop("check alpha.")}
  
  tmp_DM <- cocons::getDesignMatrix(coco.object@model.list,coco.object@data)
  
  tmp_par.pos <- tmp_DM$par.pos
  
  Hess_modified <- cocons::getModHess(coco.object = coco.object, 
                                     inv.hess = inv.hess)
  
  estims <- unlist(getEstims(coco.object)[which(unlist(lapply(tmp_par.pos,is.logical)))])[unlist(tmp_par.pos[which(unlist(lapply(tmp_par.pos, is.logical)))])]
  
  to_return <- matrix(c(estims, estims), ncol = 2, byrow = FALSE) + as.matrix(sqrt(diag(Hess_modified)), ncol = 1) %*% matrix(c(-1, 1) * stats::qnorm(alpha) ,ncol = 2)
  
  rownames(to_return) <- .cocons.updateNames(names(estims),tmp_DM$model.matrix)

  colnames(to_return) <- c(paste(( (1 - alpha)/2 )* 100, "%"),paste( (1 - ( 1 - alpha )/2) * 100, "%"))
  
  return(to_return)
}

#' Retrieves the modified inverse of the hessian
#' @description Based on the inverse of the Hessian (based on the difference parameterization for the std.dev and scale parameters), 
#' retrieves the modified inverse of the hessian (i.e. std.dev and scale).
#'
#' @usage getModHess(coco.object, inv.hess)
#' @param coco.object \code{(S4)} a fitted coco S4 object.
#' @param inv.hess \code{(matrix)} Inverse of the Hessian.
#' @returns (\code{numeric matrix}) the modified inverse of the hessian matrix
#' @author Federico Blasi
getModHess <- function(coco.object, inv.hess){
  
  Hess_modified <- inv.hess
  
  tmp_par.pos <- cocons::getDesignMatrix(coco.object@model.list,coco.object@data)$par.pos
  
  if(is.logical(tmp_par.pos$mean)){
    number_mean <- sum(tmp_par.pos$mean)
  } else{number_mean <- 0}
  
  if(is.logical(tmp_par.pos$std.dev)){
    number_std.dev <- sum(tmp_par.pos$std.dev)
  } else{number_std.dev <- 0}
  
  if(is.logical(tmp_par.pos$scale)){
    number_scale <- sum(tmp_par.pos$scale)
  } else{number_scale <- 0}
  
  if(is.logical(tmp_par.pos$std.dev) & is.logical(tmp_par.pos$scale)){
    
    to_go_over <- which(tmp_par.pos$std.dev & tmp_par.pos$scale)
    
    for(ii in 1:length(tmp_par.pos$std.dev)){
      if(tmp_par.pos$std.dev[ii] & tmp_par.pos$scale[ii]){
        
        cum_pos_std.dev <- sum(tmp_par.pos$std.dev[1:ii])
        cum_pos_scale <- sum(tmp_par.pos$scale[1:ii])
        
        Hess_modified[number_mean + cum_pos_std.dev, number_mean + cum_pos_std.dev] <- matrix(c(0.5, 0.5), ncol = 2) %*% inv.hess[c(number_mean + cum_pos_std.dev, number_mean + number_std.dev + cum_pos_scale),
                                                                                                                        c(number_mean + cum_pos_std.dev, number_mean + number_std.dev + cum_pos_scale)] %*% t(matrix(c(0.5, 0.5), ncol = 2))
        Hess_modified[number_mean + number_std.dev + cum_pos_scale, number_mean + number_std.dev + cum_pos_scale] <- matrix(c(0.5, -0.5), ncol = 2) %*% inv.hess[c(number_mean + cum_pos_std.dev, number_mean + number_std.dev + cum_pos_scale),
                                                                                                                                                       c(number_mean + cum_pos_std.dev, number_mean + number_std.dev + cum_pos_scale)] %*% t(matrix(c(0.5, -0.5), ncol = 2))
        
        
      }
    }
  }
  
  return(Hess_modified)
}

#' Retrieve the loglikelihood value
#' @description Retrieve the loglikelihood value from a fitted coco object.
#'
#' @usage getLoglik(coco.object)
#' @param coco.object \code{(S4)} a fitted coco S4 object.
#' @returns (\code{numeric}) wrap for value from a OptimParallel object 
#' @author Federico Blasi
getLoglik <- function(coco.object){
  return(coco.object@output$value)
}

#' Retrieve BIC
#' @description Retrieve BIC from a fitted coco object.
#' 
#' @usage getBIC(coco.object)
#' @param coco.object \code{(S4)} a fitted coco S4 object.
#' @returns (\code{numeric}) the associated BIC value
#' @author Federico Blasi
#' 
getBIC <- function(coco.object){
  
  if(length(coco.object@output) == 0){
    stop("object has not been fitted yet.")
  }
  
  return( coco.object@output$value +  .cocons.get.npars(coco.object) * log(dim(coco.object@data)[1]))
}

#' Retrieve AIC
#' @description Retrieve the Akaike information criterion from a fitted coco object.
#' 
#' @usage getAIC(coco.object)
#' @param coco.object \code{(S4)} a fitted coco S4 object.
#' @returns (\code{numeric}) the associated AIC value
#' @author Federico Blasi
#' 
getAIC <- function(coco.object){
  
  if(length(coco.object@output) == 0){
    stop("object has not been fitted yet.")
  }
  
  return( coco.object@output$value +  2 * log(dim(coco.object@data)[1]))
}

#' Retrieve estimates from a fitted coco object
#' @description Retrieve estimates from a fitted coco object.
#' 
#' @usage getEstims(coco.object)
#' @param coco.object \code{(S4)} a fitted coco S4 object.
#' @returns (\code{list}) a list with the estimates parameters for the different aspects
#' @author Federico Blasi
#' 
getEstims <- function(coco.object){
  
  return(cocons::getModelLists(coco.object@output$par, 
                              par.pos = getDesignMatrix(coco.object@model.list,
                                                        data = coco.object@data)$par.pos,
                              type = "diff"))
}

#' Fast and simple standardization for the design matrix.
#' @description Centers and scale the design matrix.
#' @usage getScale(x, mean.vector = NULL, sd.vector = NULL)
#' @param x \code{(S4) or (matrix)} a coco object, or a n x p matrix with covariate information to introduce, 
#' where the first column is a column of ones.
#' @param mean.vector \code{(numeric vector)} if provided, it centers covariates based on this information.
#' @param sd.vector \code{(numeric vector)} if provided, it scales covariates based on this information.
#' @returns (\code{list}) a list with a scaled design matrix of dimension n x (p+1), and a set of mean and sd vectors 
#' employed to scale the matrix
#' @author Federico Blasi
#' 
getScale <- function(x, mean.vector = NULL, sd.vector = NULL){
  
  if(methods::is(x, "coco")){
    
    x_std <- getDesignMatrix(x@model.list, x@data)$model.matrix
    
    if(is.null(mean.vector)){
      mean.vector <- apply(x_std, MARGIN = 2, base::mean)
      mean.vector[1] <- 0
    }
    
    if(is.null(sd.vector)){
      sd.vector <- apply(x_std, MARGIN = 2, stats::sd)
      sd.vector[1] <- 1
    }
    
    if(dim(x_std)[2] == 1){
      return(list("std.covs" = x_std,
                  "mean.vector" = mean.vector,
                  "sd.vector" = sd.vector))
    }
    
    for(ii in 2:dim(x_std)[2]){
      
      x_std[, ii] <- (x_std[, ii] - mean.vector[ii]) / sd.vector[ii]
      
    }
    
    return(list("std.covs" = x_std,
                "mean.vector" = mean.vector,
                "sd.vector" = sd.vector))
    
  } else{
    
    if(is.null(mean.vector)){
      mean.vector <- apply(x, MARGIN = 2, base::mean)
      mean.vector[1] <- 0
    }
    
    if(is.null(sd.vector)){
      sd.vector <- apply(x, MARGIN = 2, stats::sd)
      sd.vector[1] <- 1
    }
    
    if(dim(x)[2] == 1){
      return(list("std.covs" = x,
                  "mean.vector" = mean.vector,
                  "sd.vector" = sd.vector))
    }
    
    for(ii in 2:dim(x)[2]){
      
      x[,ii] <- (x[,ii] - mean.vector[ii]) / sd.vector[ii]
      
    }
    
    return(list("std.covs" = x,
                "mean.vector" = mean.vector,
                "sd.vector" = sd.vector))
  }
}

#' Create an efficient design matrix based on a list of aspect models
#' @description Creates a unique design matrix based on model specification for
#' each of the different potentially spatially varying aspects. 
#'
#' @usage getDesignMatrix(model.list, data)
#' @param model.list \code{(list)} a list of formulas, one for each source of nonstationarity, specifying the
#'                   models.
#' @param data \code{(data.frame)} a data.frame.
#' @return (\code{list}) a list with two elements: a design matrix of dimension 
#' (n x p), and a par.pos object, indexing columns of the design matrix to each of the spatially-varying functions.
#' @author Federico Blasi
#' 
getDesignMatrix <- function(model.list, data){
  
  first_formula <- NA_real_
  
  # search for the first formula to start the elaboration of the design matrix
  for(ii in 1:length(model.list)){
    if(is.formula(model.list[[ii]])){
      tmp_dense_formula <- model.list[[ii]]
      first_formula <- ii
      break
    }
  }
  
  # no formula detected, therefore no design matrix
  
  if(is.na(first_formula)){
    warning("No formula detected")
    return(0)
  }
  
  tmp_index_not_fixed <- unlist(lapply(model.list, is.formula))
  
  tmp_all_terms <- lapply(model.list[tmp_index_not_fixed], 
                          FUN = function(x) attr(stats::terms(x), which = "term.labels"))
  
  tmp_all_terms <- tmp_all_terms[unlist(lapply(tmp_all_terms,function(x) length(x) > 0))]
  
  # Just all non-fixed intercepts? ### here fix this
  
  if(all(unlist(lapply(tmp_all_terms, function(x) identical(x, character()))))){
    to_return_model.matrix <- rep(1, dim(data)[1])
    par.pos_list <- list()
    
    for(ii in 1:length(model.list)){
      if(is.formula(model.list[[ii]])){
        par.pos_list[[ii]] <- TRUE} else{par.pos_list[[ii]] <- model.list[[ii]]}
    }
    
    base::names(par.pos_list) <- base::names(model.list)
    
    return(list("model.matrix" = stats::model.matrix(tmp_dense_formula, data),
                "par.pos" = par.pos_list))
  }
  
  ####
  
  # If not then,
  
  tmp_dense_formula <- tmp_all_terms[[1]]
  
  if(length(tmp_all_terms) > 1){
    for(ii in 2:length(tmp_all_terms)){
      tmp_dense_formula <- union(tmp_dense_formula, tmp_all_terms[[ii]])
    }
  }
  
  tmp_intercept_yes <- any(unlist(lapply(model.list[tmp_index_not_fixed],
                                         FUN = function(x) ifelse(attr(stats::terms(x),"intercept") == 1, TRUE, FALSE))))
  
  tmp_dense_formula <- stats::as.formula(paste0(ifelse(tmp_intercept_yes, "~1+", "~0 +"),
                                                paste(tmp_dense_formula, collapse = "+")))
  
  # create design matrix
  to_return_model.matrix <- stats::model.matrix(tmp_dense_formula, data)
  
  # now indexing to create theta_list
  
  par.pos_list <- list()
  
  for(ii in 1:length(model.list)){
    
    # check if the aspect is fixed
    if(!is.formula(model.list[[ii]])){
      par.pos_list[[ii]] <- model.list[[ii]][1]
      next
    }
    
    # check if for a specific variable list is empty w/ no intercept
    if((attr(stats::terms(model.list[[ii]]), "variables") == "list()") & (attr(stats::terms(model.list[[ii]]), "intercept") == 0)) {
      par.pos_list[[ii]] <-  logical(dim(to_return_model.matrix)[2]) # return all 0's
      next
    }
    
    # check if no variables but intercept
    if((attr(stats::terms(model.list[[ii]]),"variables") == "list()") & (attr(stats::terms(model.list[[ii]]),"intercept") == 1)){
      par.pos_list[[ii]] <-  logical(dim(to_return_model.matrix)[2])
      par.pos_list[[ii]][1] <- TRUE
      next
    } 
    
    # now assuming we have variables, create the indexing
    tmp_lapply <- lapply(labels(stats::terms(model.list[[ii]])),function(x) {x == attributes(to_return_model.matrix)$dimnames[[2]]})
    par.pos_list[[ii]] <-  apply(matrix(unlist(tmp_lapply), ncol = dim(to_return_model.matrix)[2],byrow = T)
                                 , MARGIN = 2, FUN = function(x) any(x))
    
    if((base::attr(stats::terms(model.list[[ii]]), "intercept") == 1)){
      par.pos_list[[ii]][1] <- TRUE
    }
  }
  
  base::names(par.pos_list) <- base::names(model.list)
  
  return(list("model.matrix" = stats::model.matrix(tmp_dense_formula, data),
              "par.pos" = par.pos_list)
  )
}

#' Builds the necessary input for building covariance matrices
#' @description Returns a list of parameter vectors for each of the aspects.
#'
#' @usage getModelLists(theta, par.pos, type = 'diff')
#' @param theta \code{(numeric vector)} a vector of length p, where p is the number of parameters for
#' each of the models.
#' @param par.pos \code{(list)} a list detailing in which position of each aspect the elements
#' of theta should be placed. Expected to be par.pos output of \link{getDesignMatrix}.
#' @param type \code{(character)} whether parameters are related to a classical parameterization ('classic') or
#' a difference parameterization 'diff' . Default set to 'diff'.
#' @returns (\code{list}) a list of different spatial aspects and mean required for the cov.rns functions
#' @author Federico Blasi
#' 
getModelLists <- function(theta, par.pos, type = "diff"){
  
  # add a check here length theta and par.pos
  
  list_pars <- list()
  acum_pars <- 0
  
  length_logical <- max(unlist(lapply(par.pos, length)))
  
  for(ii in 1:length(par.pos)){
    if(!is.logical(par.pos[[ii]])){
      list_pars[[ii]] <- double(length_logical)
      list_pars[[ii]][1] <- par.pos[[ii]][1]
      next
    }
    
    list_pars[[ii]] <- double(length_logical)
    list_pars[[ii]][par.pos[[ii]]] <- theta[(acum_pars + 1):(acum_pars + sum(par.pos[[ii]]))]
    acum_pars <- acum_pars + sum(par.pos[[ii]])
  }
  
  names(list_pars) <- names(par.pos)
  
  if(type == "classic"){
    return(list_pars)
  }
  
  if(type == "diff"){
    
    tmp_info <- lapply(par.pos, FUN = is.logical)
    
    length_strange <- mean(unlist(lapply(par.pos[which(unlist(tmp_info))], FUN = length))) # should be only variance and range but at the end all should have the same length
    
    list_pars_temp <- list_pars
    
    if(is.logical(par.pos$std.dev) & is.logical(par.pos$scale)){
      for(ii in 1:length_strange){
        if(par.pos$std.dev[ii] & par.pos$scale[ii]){
          list_pars_temp$std.dev[ii] <- (list_pars$std.dev[ii] + list_pars$scale[ii]) / 2
          list_pars_temp$scale[ii] <- (list_pars$std.dev[ii] - list_pars$scale[ii]) / 2
        }
      }
    }
    
    return(list_pars_temp)  
  }
}

#' Simple build of boundaries
#' @description provides a set of upper and lower bounds for the L-BFGS-B routine
#' 
#' @usage getBoundariesV4(coco.object, lower.bound, upper.bound)
#' @param coco.object \code{(S4)} a coco.object.
#' @param lower.bound \code{(numeric vector)} lower bound for covariate-driven parameters of the nonstationary covariance function.
#' @param upper.bound \code{(numeric vector)} upper bound for covariate-driven parameters of the nonstationary covariance function.
#' @returns (\code{list}) a list with boundaries and simple init values for the optim L-BFGS-B routine
#' @author Federico Blasi
#' 
getBoundariesV4 <- function(coco.object, lower.bound = 2, upper.bound = 2){
  
  designMatrix <- getDesignMatrix(model.list = coco.object@model.list, data = coco.object@data)
  
  term_labels <- paste("z", "~" , paste(coco.object@model.list$mean)[2])
  
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
  
  if(is.formula(coco.object@model.list$mean)){
    coefs_lm <- stats::coef(stats::lm(term_labels,cbind.data.frame('z' = coco.object@z,mod_DM)))
  } else{
    coefs_lm <- c()
  }
  
  boundaries_B <- getBoundariesV2(coco.object = coco.object,
                                  mean.limits = c(-Inf, 0, Inf),
                                  std.dev.limits = c(-lower.bound, 0, upper.bound),
                                  scale.limits = c(-lower.bound, 0, upper.bound),
                                  aniso.limits =  c(-lower.bound, 0, upper.bound),
                                  tilt.limits =  c(-lower.bound, 0, upper.bound),
                                  smooth.limits = c(-lower.bound, 0, upper.bound),
                                  nugget.limits = c(-lower.bound, 0, upper.bound))
  
  if(is.formula(coco.object@model.list$mean)){
    boundaries_B$theta_init[1:length(coefs_lm)] <- coefs_lm
  }
  
  first_sd <- which(names(boundaries_B$theta_init) == "std.dev.limits")[1]
  n_var <- length(which(names(boundaries_B$theta_init) == "std.dev.limits")) - 1
  
  first_range <- which(names(boundaries_B$theta_init) == "scale.limits")[1]
  n_range <- length(which(names(boundaries_B$theta_init) == "scale.limits")) - 1
  
  first_smooth <- which(names(boundaries_B$theta_init) == "smooth.limits")[1]
  n_smooth <- length(which(names(boundaries_B$theta_init) == "smooth.limits")) - 1
  
  boundaries_B$theta_init[first_range] <- (log(stats::sd(coco.object@z)) - log(stats::sd(c(stats::dist(coco.object@locs)))))/2
  boundaries_B$theta_init[first_sd] <- (log(stats::sd(coco.object@z)) + log(stats::sd(c(stats::dist(coco.object@locs)))))/2
  
  boundaries_B$theta_upper[c(first_sd, first_range)] <- c(3 * abs(boundaries_B$theta_init[first_sd]), 3 * abs(boundaries_B$theta_init[first_range]))
  boundaries_B$theta_lower[c(first_sd, first_range)] <- c(-3 * abs(boundaries_B$theta_init[first_sd]), -3 * abs(boundaries_B$theta_init[first_range]))
  
  boundaries_B$theta_upper[first_smooth] <- 3
  boundaries_B$theta_lower[first_smooth] <- -3.5

  boundaries_B$theta_init[1] <- mean(coco.object@z)
  boundaries_B$theta_upper[1] <- boundaries_B$theta_init[1] + 3 * boundaries_B$theta_init[1]
  boundaries_B$theta_lower[1] <- boundaries_B$theta_init[1] - 3 * boundaries_B$theta_init[1]
  
  return(boundaries_B)
}

#' Simple build of boundaries
#' @description provides a generic set of upper and lower bounds for the L-BFGS-B routine
#' 
#' @usage getBoundaries(x, lower.value, upper.value)
#' @param x \code{(S4) or (list)} a coco.object or a par.pos list (as output from \link{getDesignMatrix})
#' @param lower.value \code{(numeric vector)} if provided, provides a vector filled with values lower.value. 
#' @param upper.value \code{(numeric vector)} if provided, provides a vector filled with values upper.value.
#' @returns (\code{list}) a list with boundaries and simple init values for the optim L-BFGS-B routine
#' @author Federico Blasi
#' 
getBoundaries <- function(x, lower.value = -2, upper.value = 2){
  
  if(upper.value < lower.value){stop("upper.value lower than lower.value")}
  
  if(methods::is(x,"coco")){
    
    tmp_par.pos <- cocons::getDesignMatrix(model.list = x@model.list,
                                          data = x@data)$par.pos
    
    tmp_index <- lapply(tmp_par.pos, FUN = is.logical)
    
    theta_init <- c(rep(0, sum(unlist(lapply(tmp_par.pos, sum))[which(tmp_index == TRUE)])))
    
    theta_upper <- rep(upper.value, length(theta_init))
    theta_lower <- rep(lower.value, length(theta_init))
    
    return(list("theta_init" = theta_init,
                "theta_upper" = theta_upper,
                "theta_lower" = theta_lower))
    
  } else {stop("provide an coco object.")}
  
}

#' Simple build of boundaries (v2)
#' @description provides a generic set of upper and lower bounds for the L-BFGS-B routine
#'
#' @usage getBoundariesV2(coco.object, mean.limits, std.dev.limits, 
#' scale.limits, aniso.limits, tilt.limits, smooth.limits, nugget.limits)
#' @param coco.object \code{(S4)} a coco object.
#' @param mean.limits \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param std.dev.limits \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param scale.limits \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param aniso.limits \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param tilt.limits \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param smooth.limits \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param nugget.limits \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @returns (\code{list}) a list with boundaries for the optim L-BFGS-B routine
#' @author Federico Blasi
#' 
getBoundariesV2 <- function(coco.object, 
                            mean.limits,
                            std.dev.limits,
                            scale.limits,
                            aniso.limits,
                            tilt.limits,
                            smooth.limits,
                            nugget.limits){
  
  if(!is.null(coco.object)){
    if(methods::is(coco.object,"coco")){
      
      limits <- rbind(mean.limits,
                      std.dev.limits,
                      scale.limits,
                      aniso.limits,
                      tilt.limits,
                      smooth.limits,
                      nugget.limits)
      
      tmp_par.pos <- getDesignMatrix(model.list = coco.object@model.list,
                                     data = coco.object@data)$par.pos
      
      tmp_index <- lapply(tmp_par.pos, FUN = is.logical)
      
      theta_init <- theta_upper <- theta_lower <- numeric(0)
      
      for(ii in 1:length(tmp_index)){
        
        if(!tmp_index[[ii]]){next}
        
        theta_lower <- append(theta_lower,rep(limits[ii,1],sum(tmp_par.pos[[ii]])))
        theta_init <- append(theta_init,rep(limits[ii,2],sum(tmp_par.pos[[ii]])))
        theta_upper <- append(theta_upper,rep(limits[ii,3],sum(tmp_par.pos[[ii]])))
        
      }
      
      return(list("theta_init" = theta_init,
                  "theta_upper" = theta_upper,
                  "theta_lower" = theta_lower))
      
    }
    
  } else{stop("coco.object required")}
}

#' Simple build of boundaries (v3)
#' @description provides a generic set of upper and lower bounds for the L-BFGS-B routine
#'
#' @usage getBoundariesV3(coco.object, mean.limits, global.lower, 
#' std.dev.max.effects, 
#' scale.max.effects, aniso.max.effects, tilt.max.effects, 
#' smooth.max.effects, nugget.max.effects)
#' @param coco.object \code{(S4)} a coco object.
#' @param mean.limits \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param global.lower \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param std.dev.max.effects \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param scale.max.effects \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param aniso.max.effects \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param tilt.max.effects \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param smooth.max.effects \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @param nugget.max.effects \code{(numeric vector)} a vector of c(lower,init,upper) values for the associated param.
#' @returns (\code{list}) a list with boundaries for the optim L-BFGS-B routine
#' @author Federico Blasi
#' 
getBoundariesV3 <- function(coco.object,
                            mean.limits,
                            global.lower,
                            std.dev.max.effects,
                            scale.max.effects,
                            aniso.max.effects,
                            tilt.max.effects,
                            smooth.max.effects,
                            nugget.max.effects){
  
  if(!is.null(coco.object)){
    if(methods::is(coco.object,"coco")){
      
      # lower bounds for global effects
      
      emp_var <- stats::var(coco.object@z)
      
      max_dist <- max(as.matrix(stats::dist(coco.object@locs))) * 0.7
      
      limits <- rbind(mean.limits,
                      std.dev.max.effects,
                      scale.max.effects,
                      aniso.max.effects,
                      tilt.max.effects,
                      smooth.max.effects,
                      nugget.max.effects)
      
      tmp_par.pos <- getDesignMatrix(model.list = coco.object@model.list,
                                     data = coco.object@data)$par.pos
      
      tmp_index <- lapply(tmp_par.pos, FUN = is.logical)
      
      theta_init <- theta_upper <- theta_lower <- numeric(0)
      
      # trend
      
      theta_lower <- append(theta_lower, rep(-mean.limits,sum(tmp_par.pos[[1]])))
      theta_init <- append(theta_init, rep(0,sum(tmp_par.pos[[1]])))
      theta_upper <- append(theta_upper, rep(mean.limits, sum(tmp_par.pos[[1]])))
      
      
      for(ii in c(2:4)){
        
        if(!tmp_index[[ii]]){next}
        
        theta_lower <- append(theta_lower, rep(ifelse(limits[ii,1] > 1 ,
                                                      yes = log(1/limits[ii,1]),
                                                      no = log(limits[ii,1])), sum(tmp_par.pos[[ii]])))
        theta_init <- append(theta_init, rep(0,sum(tmp_par.pos[[ii]])))
        theta_upper <- append(theta_upper, rep(ifelse(limits[ii, 1] > 1 ,
                                                      yes = log(limits[ii, 1]),
                                                      no = log(1/limits[ii, 1])), sum(tmp_par.pos[[ii]])))
        
      }
      
      # tilt
      
      if(is.logical(tmp_par.pos$tilt)){
        
        names(tilt.max.effects) <- "tilt.max.effects"
        theta_lower <- append(theta_lower, rep(-tilt.max.effects,sum(tmp_par.pos[[5]])))
        theta_init <- append(theta_init, rep(0,sum(tmp_par.pos[[5]])))
        theta_upper <- append(theta_upper, rep(tilt.max.effects, sum(tmp_par.pos[[5]])))
        
      }
      
      # smtns
      
      if(is.logical(tmp_par.pos$smooth)){
        
        names(smooth.max.effects) <- "smooth.max.effects"
        theta_lower <- append(theta_lower, rep(-smooth.max.effects, sum(tmp_par.pos[[6]])))
        theta_init <- append(theta_init, rep(0,sum(tmp_par.pos[[6]])))
        theta_upper <- append(theta_upper, rep(smooth.max.effects, sum(tmp_par.pos[[6]])))
        
      }
      
      # nugget
      
      if(is.logical(tmp_par.pos$nugget)){
        
        names(nugget.max.effects) <- "nugget.max.effects"
        theta_lower <- append(theta_lower, rep(log(1/nugget.max.effects), sum(tmp_par.pos[[7]])))
        theta_init <- append(theta_init, rep(0,sum(tmp_par.pos[[7]])))
        theta_upper <- append(theta_upper, rep(log(nugget.max.effects), sum(tmp_par.pos[[7]])))
        
      }
      
      first_var <- which(names(theta_lower) == "std.dev.max.effects")[1]
      first_range <- which(names(theta_lower) == "scale.max.effects")[1]
      first_smooth <- which(names(theta_lower) == "smooth.max.effects")[1]
      
      theta_lower[first_var] <- 2 * log(global.lower)
      theta_lower[first_range] <- log(global.lower) - log(max_dist)
      
      theta_upper[first_var] <- log(3 * emp_var) + log(max_dist)
      theta_upper[first_range] <- log(3 * emp_var) - log(global.lower)
      
      theta_init[first_var] <- mean(c(theta_upper[first_var],theta_lower[first_var]))
      theta_init[first_range] <- mean(c(theta_upper[first_range],theta_lower[first_range]))
      
      theta_init[first_smooth] <- theta_lower[first_smooth]
      
      return(list("theta_init" = theta_init,
                  "theta_upper" = theta_upper,
                  "theta_lower" = theta_lower))
      
    }
    
  } else{stop("coco object required")}
}

#' getHessian
#' @description numerically approximate the Hessian. Hessians of parameters based on "pml" are based on full likelihoods.
#' @usage getHessian(coco.object, ncores = parallel::detectCores() - 1, 
#' eps = .Machine$double.eps^(1/4))
#' @param coco.object \code{(S4)} a fitted coco object.
#' @param ncores \code{(integer)} number of cores used for the computation.
#' @param eps \code{(numeric)} ...
#' @returns (\code{numeric matrix}) a symmetric matrix pxp of the approximated (observed) Hessian
#' @author Federico Blasi
getHessian <- function(coco.object, ncores = parallel::detectCores() - 1, 
                       eps = .Machine$double.eps^(1/4)){
  
  if(coco.object@type == "dense"){
    
    if(coco.object@info$optim.type == "reml"){stop("reml hessian not implemented yet.")}
    
    p <- base::length(coco.object@output$par)
    H <- base::matrix(NA, p, p)
    
    index_positions <- base::matrix(ncol = 2)
    
    for(jj in 1:(p)){
      for(ii in (jj):p){
        index_positions <- base::rbind(index_positions, c(jj, ii))
      }
    }
    
    index_positions <- index_positions[-1, ]
    
    Design_Matrix <- cocons::getDesignMatrix(model.list = coco.object@model.list, 
                                            data = coco.object@data)
    par.pos <- Design_Matrix$par.pos
    x_covariates <- Design_Matrix$model.matrix
    
    coco_x_std <- cocons::getScale(x_covariates, 
                                   mean.vector = coco.object@info$mean.vector,
                                   sd.vector = coco.object@info$sd.vector)
    
    n <- dim(coco.object@z)[1]
    
    x_covariates <- coco_x_std$std.covs
    
    if(coco.object@info$optim.type == "pml"){
      
      f00 <- cocons::GetNeg2loglikelihood(coco.object@output$par, par.pos = par.pos,
                                   locs = coco.object@locs,
                                   x_covariates = x_covariates, 
                                   smooth.limits = coco.object@info$smooth.limits, 
                                   n = n, 
                                   z = coco.object@z,
                                   lambda = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg))
      
    } else{f00 <- coco.object@output$value}
    
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::setDefaultCluster(cl = cl)
    parallel::clusterEvalQ(cl, library("cocons"))
    
    parallel::clusterExport(cl = cl, list("x_covariates", "coco.object",
                                          "n", "par.pos", "eps", "f00"),
                            envir = environment())
    
    vector_responses <- parallel::parApply(cl = cl, 
                                           index_positions, 
                                           MARGIN = 1, 
                                           FUN = function(x){
                                             
                                             t01 <- t10 <- t11 <- coco.object@output$par

                                             t01[x[1]] <- t01[x[1]] + eps
                                             t10[x[2]] <- t10[x[2]] + eps
                                             t11[x[1]] <- t11[x[1]] + eps
                                             t11[x[2]] <- t11[x[2]] + eps
                                             
                                             f01 <- cocons::GetNeg2loglikelihood(t01, par.pos = par.pos,
                                                                                    locs = coco.object@locs,
                                                                                    x_covariates = x_covariates, 
                                                                                    smooth.limits = coco.object@info$smooth.limits, 
                                                                                    n = n, 
                                                                                    z = coco.object@z,
                                                                                    lambda = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg))
                                             
                                             f10 <- cocons::GetNeg2loglikelihood(t10, par.pos = par.pos,
                                                                                    locs = coco.object@locs,
                                                                                    x_covariates = x_covariates, 
                                                                                    smooth.limits = coco.object@info$smooth.limits, 
                                                                                    n = n, 
                                                                                    z = coco.object@z,
                                                                                    lambda = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg))
                                             
                                             f11 <- cocons::GetNeg2loglikelihood(t11, par.pos = par.pos,
                                                                                    locs = coco.object@locs,
                                                                                    x_covariates = x_covariates, 
                                                                                    smooth.limits = coco.object@info$smooth.limits, 
                                                                                    n = n, 
                                                                                    z = coco.object@z, 
                                                                                    lambda = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg))
                                             
                                             0.5 * ((f11 - f01 - f10 + f00) / (eps*eps))
                                           })
    
    count_index <- 0
    
    for(jj in 1:p){
      for(ii in jj:p){
        
        count_index <- count_index + 1
        H[jj,ii] <- vector_responses[count_index]
        
      }
    }
    
    H[is.na(H)] <- 0
    H <- H + t(H)
    diag(H) <- diag(H)/2
    
    return(H)
  }
  
  if(coco.object@type == "sparse"){
    
    p <- base::length(coco.object@output$par)
    H <- base::matrix(NA, p, p)
    
    index_positions <- base::matrix(ncol = 2)
    
    for(jj in 1:(p)){
      for(ii in (jj):p){
        index_positions <- base::rbind(index_positions, c(jj, ii))
      }
    }
    
    index_positions <- index_positions[-1, ]
    
    Design_Matrix <- cocons::getDesignMatrix(model.list = coco.object@model.list, 
                                            data = coco.object@data)
    par.pos <- Design_Matrix$par.pos
    x_covariates <- Design_Matrix$model.matrix
    
    coco_x_std <- cocons::getScale(x_covariates, 
                                   mean.vector = coco.object@info$mean.vector,
                                   sd.vector = coco.object@info$sd.vector)
    
    n <- dim(coco.object@z)[1]
    x_covariates <- coco_x_std$std.covs
    
    ref_taper <- coco.object@info$taper(
      spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL),
      theta = c(coco.object@info$delta, 1)
    )
    
    cholS <- spam::chol.spam(ref_taper)
    
    if(coco.object@info$optim.type == "pml"){
      
      f00 <- cocons::GetNeg2loglikelihoodTaper(theta = coco.object@output$par,
                                               par.pos = par.pos,
                                               ref_taper = ref_taper,
                                               locs = coco.object@locs,
                                               x_covariates = x_covariates,
                                               smooth.limits = coco.object@info$smooth.limits,
                                               cholS = cholS,
                                               z = coco.object@z,
                                               n = n,
                                               lambda = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg))
      
    } else{f00 <- coco.object@output$value}

    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::setDefaultCluster(cl = cl)
    parallel::clusterEvalQ(cl, library("cocons"))
    
    parallel::clusterExport(cl = cl, list("x_covariates", "coco.object",
                                          "n", "par.pos", "eps",
                                          "ref_taper","cholS","f00"),
                            envir = environment())
    
    vector_responses <- parallel::parApply(cl = cl, 
                                           index_positions, 
                                           MARGIN = 1, 
                                           FUN = function(x){
                                             
                                             t11 <- t10 <- t01 <- coco.object@output$par

                                             t01[x[1]] <- t01[x[1]] + eps
                                             t10[x[2]] <- t10[x[2]] + eps
                                             t11[x[1]] <- t11[x[1]] + eps
                                             t11[x[2]] <- t11[x[2]] + eps
                                             
                                             f01 <- cocons::GetNeg2loglikelihoodTaper(theta = t01,
                                                                                          par.pos = par.pos,
                                                                                          ref_taper = ref_taper,
                                                                                          locs = coco.object@locs,
                                                                                          x_covariates = x_covariates,
                                                                                          smooth.limits = coco.object@info$smooth.limits,
                                                                                          cholS = cholS,
                                                                                          z = coco.object@z,
                                                                                          n = n,
                                                                                          lambda = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg))
                                             
                                             f10 <- cocons::GetNeg2loglikelihoodTaper(theta = t10,
                                                                                          par.pos = par.pos,
                                                                                          ref_taper = ref_taper,
                                                                                          locs = coco.object@locs,
                                                                                          x_covariates = x_covariates,
                                                                                          smooth.limits = coco.object@info$smooth.limits,
                                                                                          cholS = cholS,
                                                                                          z = coco.object@z,
                                                                                          n = n,
                                                                                          lambda = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg))
                                             
                                             f11 <- cocons::GetNeg2loglikelihoodTaper(theta = t11,
                                                                                          par.pos = par.pos,
                                                                                          ref_taper = ref_taper,
                                                                                          locs = coco.object@locs,
                                                                                          x_covariates = x_covariates,
                                                                                          smooth.limits = coco.object@info$smooth.limits,
                                                                                          cholS = cholS,
                                                                                          z = coco.object@z,
                                                                                          n = n,
                                                                                          lambda = c(coco.object@info$lambda.Sigma, coco.object@info$lambda.betas, coco.object@info$lambda.reg))
                                             
                                             0.5 * ((f11 - f01 - f10 + f00) / (eps*eps))
                                             
                                           })
    
    count_index <- 0
    
    for(jj in 1:p){
      for(ii in jj:p){
        
        count_index <- count_index + 1
        H[jj,ii] <- vector_responses[count_index]
        
      }
    }
    
    H[is.na(H)] <- 0
    H <- H + t(H)
    diag(H) <- diag(H)/2
    
    return(H)
    
  }
 
  
}
