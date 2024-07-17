# New ones

#' Covariance matrix from a fitted svcov object
#' @description Retrieves the associated covariance matrix from a fitted svcov object
#'
#' @usage getCovmatrix(svcov.object, type = 'global', index = NULL)
#' @param svov.object an svcov class (fitted) object
#' @param type whether 'global' to retrieve the regular covariance matrix, or 'local' to retrieve global covariance 
#' based on the local aspect of a specific location (not implemented yet)
#' @param index index to perform local covariance matrix (not implemented yet)
#' @returns a vector with the adjusted trend
#' @author Federico Blasi
getCovMatrix <- function(svcov.object, type = 'global', index = NULL){
  
  x_covs <- svcov::getScale(svcov.object)$std.covs
  
  par.pos <- getDesignMatrix(svcov.object@model.list, data = svcov.object@data)$par.pos
  
  if(svcov.object@type == 'dense'){
    
    if(type == 'global'){
      theta_list <- svcov::getModelLists(svcov.object@output$par,par.pos = par.pos, 
                                         type = 'diff')
      
      return(svcov::cov_rns(theta = theta_list,
                            locs = svcov.object@locs,
                            x_covariates = x_covs,
                            smooth_limits = svcov.object@info$smooth_limits))
    }
    
    if(type == 'local'){
      
    }
    
  }
  
  if(svcov.object@type == 'sparse'){
    
    if(type == 'global'){
      
      
      theta_list <- svcov::getModelLists(svcov.object@output$par,par.pos = par.pos, 
                                         type = 'diff')
      
      # taper
      ref_taper <- svcov.object@info$taper(
        spam::nearest.dist(svcov.object@locs, delta = svcov.object@info$delta, upper = NULL),
        theta = c(svcov.object@info$delta, 1)
      )
      
      ref_taper@entries <- ref_taper@entries * cov_rns_taper_optimized_range(theta = theta_list[-1], 
                                                                                 locs = svcov.object@locs, 
                                                                                 x_covariates =  x_covs, 
                                                                                 colindices = ref_taper@colindices, 
                                                                                 rowpointers = ref_taper@rowpointers,
                                                                                 smooth_limits =  svcov.object@info$smooth_limits)
      
      return(as.matrix(ref_taper))
      
    }
    
    if(type == 'local'){
      
      
    }
    
  }
}

#' Based on a set of predictions retrieves the Logrank
#' @description Retrieves the estimated spatial effects of the spatial structure
#'
#' @usage getLogrank(svcov.object, type = 'global', index = NULL)
#' @param svov.object an svcov class (fitted) object
#' @param type whether 'global' to retrieve the regular covariance matrix, or 'local' to retrieve global covariance 
#' based on the local aspect of a specific location (not implemented yet)
#' @param index index to perform local covariance matrix (not implemented yet)
#' @returns a vector with the adjusted trend
#' @author Federico Blasi
getLogRank <- function(z.pred, mean.pred, sd.pred){
  
  return( (log(2 * pi) + ((z.pred - mean.pred) / sd.pred)^2)/2 + log(sd.pred))
  
}

#' Based on a set of predictions retrieves the Logrank
#' @description Retrieves the estimated spatial effects of the spatial structure
#'
#' @usage getCRPS(svcov.object, type = 'global', index = NULL)
#' @param svov.object an svcov class (fitted) object
#' @param type whether 'global' to retrieve the regular covariance matrix, or 'local' to retrieve global covariance 
#' based on the local aspect of a specific location (not implemented yet)
#' @param index index to perform local covariance matrix (not implemented yet)
#' @returns a vector with the adjusted trend
#' @author Federico Blasi
getCRPS <- function(z.pred, mean.pred, sd.pred){
  
  vector_tmp_z <- (mean.pred - z.pred) / sd.pred
  
  return( sd.pred * (vector_tmp_z * (2 * stats::pnorm(vector_tmp_z) - 1) + 
                                               2 * stats::dnorm(vector_tmp_z) - 1 / base::sqrt(pi)))

}

#' Retrieves the estimated spatial effects of the spatial structure
#' @description Retrieves the estimated spatial effects of the spatial structure
#'
#' @usage getSpateffects(svcov.object, type = 'global', index = NULL)
#' @param svov.object an svcov class (fitted) object
#' @param type whether 'global' to retrieve the regular covariance matrix, or 'local' to retrieve global covariance 
#' based on the local aspect of a specific location (not implemented yet)
#' @param index index to perfrom local covariance matrix (not implemented yet)
#' @returns a vector with the adjusted trend
#' @author Federico Blasi
getSpatEffects <- function(svcov.object){
  
  tmp_info <- svcov::getDesignMatrix(model.list = svcov.object@model.list, 
                                     data = svcov.object@data)
  
  theta_list <- svcov::getModelLists(
    theta = svcov.object@output$par, par.pos = tmp_info$par.pos,
    type = "diff"
  )
  
  X_std <- svcov::getScale(tmp_info$model.matrix,
                           mean.vector = svcov.object@info$mean.vector,
                           sd.vector = svcov.object@info$sd.vector
  )
  
  tp_se <- exp(0.5 * X_std$std.covs %*% theta_list$std.dev)
  tp_ga <- exp(X_std$std.covs %*% theta_list$aniso)
  tp_tl <- pi / (1 + exp(-X_std$std.covs %*% theta_list$tilt))
  tp_smooth <- (svcov.object@info$smooth_limits[2] - svcov.object@info$smooth_limits[1]) / (1 + exp(-X_std$std.covs %*% theta_list$smooth)) + svcov.object@info$smooth_limits[1]
  tp_ng <- exp(X_std$std.covs %*% theta_list$nugget)
  tp_mr_x <- sin(tp_tl) * exp(X_std$std.covs %*% theta_list$scale)
  tp_mr_y <- sin(tp_tl) * exp(X_std$std.covs %*% theta_list$scale) * exp(X_std$std.covs %*% theta_list$aniso)
  
  return(list('sd' = tp_se,
              'scale_x' = tp_mr_x,
              'scale_y' = tp_mr_y,
              'aniso' = tp_ga,
              'tilt' = cos(tp_tl),
              'smooth' = tp_smooth,
              'nugget' = tp_ng))
  
}

# getLogliktest <- function(svcov.full, svcov.reduced, alpha){
#  svcov.reduced@output$value - svcov.full@output$value
#}

#' Computes the condition number of the associated correlation matrix of the fitted svcov object
#' @description Compute the trend of the (fitted) svcov object
#'
#' @usage getCondNumber(svcov.object)
#' @param svov.object an svcov class (fitted) object
#' @returns the condition number
#' @author Federico Blasi
getCondNumber <- function(svcov.object){
  corr_mat <- cov2cor(getCovMatrix(svcov.object))
  eigen(corr_mat)$values[1] / eigen(corr_mat)$values[dim(corr_mat)[1]]
}

#' Computes the trend of the svcov object
#' @description Compute the trend of the (fitted) svcov object
#'
#' @usage getTrend(svcov.object)
#' @param svov.object an svcov class (fitted) object
#' @returns a vector with the adjusted trend
#' @author Federico Blasi
getTrend <- function(svcov.object){
  tmp_scaled <- getScale(svcov.object)$std.covs
  return(tmp_scaled %*% getEstims(svcov.object)$mean)
}

#' Compute Confidence Intervals for an svcov object
#' @description Compute confidence intervals for a (fitted) svcov object
#'
#' @usage getCIs(svcov.object, inv.hess, alpha = 0.05)
#' @param svov.object an svcov class (fitted) object
#' @param Hess the hessian 
#' @param alpha confidence level
#' @returns a matrix with lower and upper bound values
#' @author Federico Blasi
getCIs <- function(svcov.object, inv.hess, alpha = 0.05){
  
  tmp_par.pos <- svcov::getDesignMatrix(svcov.object@model.list,svcov.object@data)$par.pos
  
  Hess_modified <- svcov::getModHess(svcov.object = svcov.object, 
                                     inv.hess = inv.hess)
  
  estims <- unlist(getEstims(svcov.object)[which(unlist(lapply(tmp_par.pos,is.logical)))])[unlist(tmp_par.pos[which(unlist(lapply(tmp_par.pos, is.logical)))])]
  
  to_return <- matrix(c(estims, estims), ncol = 2, byrow = FALSE) + as.matrix(sqrt(diag(Hess_modified)), ncol = 1) %*% matrix(c(-1, 1) * stats::qnorm(1 - alpha/2) ,ncol = 2)
  
  rownames(to_return) <- names(estims)
  
  colnames(to_return) <- c(paste((alpha/2 )* 100, "%"),paste((1 - alpha/2 ) * 100, "%"))
  
  return(to_return)
}

#' Retrieves the modified inverse of the hessian
#' @description Compute confidence intervals for a (fitted) svcov object
#'
#' @usage getModHess(svcov.object, inv.hess, alpha = 0.05)
#' @param svov.object an svcov class (fitted) object
#' @param inv.hess Inverse of the Hessian 
#' @param alpha confidence level
#' @returns a matrix with lower and upper bound values
#' @author Federico Blasi
getModHess <- function(svcov.object, inv.hess){
  
  Hess_modified <- inv.hess
  
  tmp_par.pos <- svcov::getDesignMatrix(svcov.object@model.list,svcov.object@data)$par.pos
  
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
#' @description Retrieve the loglikelihood value from a fitted svcov object.
#'
#' @usage getLoglik(svcov.object)
#' @param svov.object an svcov class (fitted) object
#' @returns loglikelihood value of the OptimParallel function
#' @author Federico Blasi
getLoglik <- function(svcov.object){
  return(svcov.object@output$value)
}

#' Retrieve BIC
#' @description Retrieve BIC from a fitted svcov object
#' 
#' @usage getBIC(svcov.object)
#' @param svcov.object a fitted svcov object.
#' @returns a list with the associated BIC value
#' @author Federico Blasi
#' 
getBIC <- function(svcov.object){
  
  if(length(svcov.object@output) == 0){
    stop('object has not been fitted yet.')
  }
  
  temp_par.pos <- svcov::getDesignMatrix(svcov.object@model.list,svcov.object@data)$par.pos
  tmp_index <- lapply(temp_par.pos, FUN = is.logical)
  n_par <- sum(unlist(lapply(temp_par.pos, sum))[which(tmp_index == TRUE)])
  return( svcov.object@output$value +  n_par * log(dim(svcov.object@data)[1]))
}

#' Retrieve AIC
#' @description Retrieve the Akaike information criterion from a fitted svcov object
#' 
#' @usage getAIC(svcov.object)
#' @param svcov.object a fitted svcov object.
#' @returns a list with the associated AIC value
#' @author Federico Blasi
#' 
getAIC <- function(svcov.object){
  
  if(length(svcov.object@output) == 0){
    stop('object has not been fitted yet.')
  }
  
  temp_par.pos <- svcov::getDesignMatrix(svcov.object@model.list,svcov.object@data)$par.pos
  tmp_index <- lapply(temp_par.pos, FUN = is.logical)
  n_par <- sum(unlist(lapply(temp_par.pos, sum))[which(tmp_index == TRUE)])
  return( svcov.object@output$value +  2 * log(dim(svcov.object@data)[1]))
}

#' Retrieve estimates from a fitted svcov object
#' @description Retrieve estimates from a fitted svcov object
#' 
#' @usage getEstims(svcov.object)
#' @param svcov.object a fitted svcov object.
#' @returns a list with the estimates parameters for the different aspects
#' @author Federico Blasi
#' 
getEstims <- function(svcov.object){
  
  return(svcov::getModelLists(svcov.object@output$par, 
                              par.pos = getDesignMatrix(svcov.object@model.list,
                                                        data = svcov.object@data)$par.pos,
                              type = 'diff'))
}

#' Fast and simple standardization for the design matrix
#' @description Centers and scale the design matrix
#' @usage getScale(x, mean.vector = NULL, sd.vector = NULL)
#' @param x a svcov object, or a n x p matrix with covariate information to introduce, 
#' where the first column is a column of ones.
#' @param mean.vector if provided, it centers covariates based on this information
#' @param sd.vector if provided, it scales covariates based on this information
#' @returns a list with a scaled design matrix of dimension n x (p+1), and a set of mean and sd vectors 
#' employed to scale the matrix
#' @author Federico Blasi
#' 
getScale <- function(x, mean.vector = NULL, sd.vector = NULL){
  
  if('svcov' %in% class(x)){
    
    x_std <- getDesignMatrix(x@model.list, x@data)$model.matrix
    
    if(is.null(mean.vector)){
      mean.vector <- apply(x_std, MARGIN = 2, base::mean)
    }
    
    if(is.null(sd.vector)){
      sd.vector <- apply(x_std, MARGIN = 2, stats::sd)
      sd.vector[1] <- 1
    }
    
    if(dim(x_std)[2] == 1){
      return(list('std.covs' = x_std,
                  'mean.vector' = mean.vector,
                  'sd.vector' = sd.vector))
    }
    
    for(ii in 2:dim(x_std)[2]){
      
      x_std[, ii] <- (x_std[, ii] - mean.vector[ii]) / sd.vector[ii]
      
    }
    
    return(list('std.covs' = x_std,
                'mean.vector' = mean.vector,
                'sd.vector' = sd.vector))
    
  } else{
    
    if(is.null(mean.vector)){
      mean.vector <- apply(x, MARGIN = 2, base::mean)
    }
    
    if(is.null(sd.vector)){
      sd.vector <- apply(x, MARGIN = 2, stats::sd)
      sd.vector[1] <- 1
    }
    
    if(dim(x)[2] == 1){
      return(list('std.covs' = x,
                  'mean.vector' = mean.vector,
                  'sd.vector' = sd.vector))
    }
    
    for(ii in 2:dim(x)[2]){
      
      x[,ii] <- (x[,ii] - mean.vector[ii]) / sd.vector[ii]
      
    }
    
    return(list('std.covs' = x,
                'mean.vector' = mean.vector,
                'sd.vector' = sd.vector))
  }
}

#' Create an efficient design matrix based on a list of aspect models
#' @description Creates a unique design matrix based on model specification for
#' each of the different potentially spatially varying aspects. 
#'
#' @param model.list a list of formulas, one for each aspect, specifying the
#'                   models.
#' @param data a data.frame
#' @return a list() with two elements: a design matrix of dimension 
#' (n x p), and a par.pos object, indexing columns of the design matrix refered
#' to each aspect models.
#'
#' @examples 
#' model.list <- list(   "mean" = 0,
#'                       "std.dev" = as.formula(" ~ 1 + lati_s * long_s"),
#'                       "scale" = as.formula(" ~ 1 + elev_s"),
#'                       "aniso" = as.formula(" ~ 1 + elev_s"),
#'                       "tilt" = as.formula(" ~ 1 + elev_s"),
#'                       "smooth" = as.formula(" ~ 1"),
#'                       "nugget" = -Inf)
#'                       
#' @author Federico Blasi
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
    warning('No formula detected')
    return(0)
  }
  
  tmp_index_not_fixed <- unlist(lapply(model.list, is.formula))
  
  tmp_all_terms <- lapply(model.list[tmp_index_not_fixed], 
                          FUN = function(x) attr(stats::terms(x), which = "term.labels"))
  
  # added on 03/01/24
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
    
    return(list('model.matrix' = stats::model.matrix(tmp_dense_formula, data),
                'par.pos' = par.pos_list))
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
  
  return(list('model.matrix' = stats::model.matrix(tmp_dense_formula, data),
              'par.pos' = par.pos_list)
  )
}

#' Builds the necessary input for building covariance matrices
#' @description Returns a list of parameter vectors for each of the aspects.
#'
#' @usage getModelLists(optim$pars, par.pos)
#' @param theta a vector of length p, where p is the number of parameters for
#' each of the models
#' @param par.pos a list detailing in which position of each aspect the elements
#' of theta should be placed. Expected to be output of getDesignMatrix
#' @param type whether parameters are related to a classical parameterization ('classic') or
#' a difference parameterization 'diff' . Default set to 'diff'.
#' @returns a list() of different spatial aspects and mean required for the cov.rns functions
#' @examples 
#' model.list <- list(   "mean" = 0,
#'                       "std.dev" = as.formula(" ~ 1 + lati_s * long_s"),
#'                       "scale" = as.formula(" ~ 1 + elev_s"),
#'                       "aniso" = as.formula(" ~ 1 + elev_s"),
#'                       "tilt" = as.formula(" ~ 1 + elev_s"),
#'                       "smooth" = as.formula(" ~ 1"),
#'                       "nugget" = -Inf)
#'                       
#' @author Federico Blasi
getModelLists <- function(theta, par.pos, type = 'diff'){
  
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
  
  if(type == 'classic'){
    return(list_pars)
  }
  
  if(type == 'diff'){
    
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
#' @description provides a generic set of upper and lower bounds for the L-BFGS-B routine
#' 
#' @usage getBoundaries(par.pos, lower.value, upper.value)
#' @param x an svcov.object or a par.pos list (as output from getDesignMatrix)
#' @param lower.value if provided, provides a vector filled with values lower.value. 
#' @param upper.value if provided, provides a vector filled with values upper.value.
#' @returns a list with boundaries and simple init values for the optim L-BFGS-B routine
#' @examples 
#' model.list <- list("std.dev" = as.formula(" ~ 1 + lati_s * long_s"),
#' "scale" = as.formula(" ~ 1 + elev_s"),
#' "aniso" = as.formula(" ~ 1 + elev_s"),
#' "tilt" = as.formula(" ~ 1 + elev_s"),
#' "smooth" = as.formula(" ~ 1"),
#' "nugget" = -Inf)
#' @author Federico Blasi
#' 
getBoundaries <- function(x, lower.value, upper.value){
  
  if(upper.value < lower.value){stop('upper.value lower than lower.value')}
  
  if(class(x) == 'svcov'){
    
    tmp_par.pos <- svcov::getDesignMatrix(model.list = x@model.list,
                                          data = x@data)$par.pos
    
    tmp_index <- lapply(tmp_par.pos, FUN = is.logical)
    
    theta_init <- c(rep(0, sum(unlist(lapply(tmp_par.pos, sum))[which(tmp_index == TRUE)])))
    
    theta_upper <- rep(upper.value, length(theta_init))
    theta_lower <- rep(lower.value, length(theta_init))
    
    return(list('theta_init' = theta_init,
                'theta_upper' = theta_upper,
                'theta_lower' = theta_lower))
    
  } else {stop('provide an svcov object.')}
  
}

#' Simple build of boundaries
#' @description provides a generic set of upper and lower bounds for the L-BFGS-B routine
#'
#' @usage getBoundariesV2(par.pos, lower.value, upper.value)
#' @param svcov.object a par.pos object output from getDesignMatrix
#' @param param.limits a vector of c(lower,init,upper) values for the associated param.
#' @returns a list with boundaries for the optim L-BFGS-B routine
#' @examples 
#' model.list <- list("std.dev" = as.formula(" ~ 1 + lati_s * long_s"),
#' "scale" = as.formula(" ~ 1 + elev_s"),
#' "aniso" = as.formula(" ~ 1 + elev_s"),
#' "tilt" = as.formula(" ~ 1 + elev_s"),
#' "smooth" = as.formula(" ~ 1"),
#' "nugget" = -Inf)
#' @author Federico Blasi
#' 
getBoundariesV2 <- function(svcov.object, 
                            mean.limits,
                            std.dev.limits,
                            scale.limits,
                            aniso.limits,
                            tilt.limits,
                            smooth.limits,
                            nugget.limits){
  
  if(!is.null(svcov.object)){
    if(class(svcov.object) == 'svcov'){
      
      limits <- rbind(mean.limits,
                      std.dev.limits,
                      scale.limits,
                      aniso.limits,
                      tilt.limits,
                      smooth.limits,
                      nugget.limits)
      
      tmp_par.pos <- getDesignMatrix(model.list = svcov.object@model.list,
                                     data = svcov.object@data)$par.pos
      
      tmp_index <- lapply(tmp_par.pos, FUN = is.logical)
      
      theta_init <- theta_upper <- theta_lower <- numeric(0)
      
      for(ii in 1:length(tmp_index)){
        
        if(!tmp_index[[ii]]){next}
        
        theta_lower <- append(theta_lower,rep(limits[ii,1],sum(tmp_par.pos[[ii]])))
        theta_init <- append(theta_init,rep(limits[ii,2],sum(tmp_par.pos[[ii]])))
        theta_upper <- append(theta_upper,rep(limits[ii,3],sum(tmp_par.pos[[ii]])))
        
      }
      
      return(list('theta_init' = theta_init,
                  'theta_upper' = theta_upper,
                  'theta_lower' = theta_lower))
      
    }
    
  } else{stop('svcov.object required')}
}

getBoundariesV3 <- function(svcov.object,
                            mean.limits,
                            global.lower,
                            std.dev.max.effects,
                            scale.max.effects,
                            aniso.max.effects,
                            tilt.max.effects,
                            smooth.max.effects,
                            nugget.max.effects){
  
  if(!is.null(svcov.object)){
    if(class(svcov.object) == 'svcov'){
      
      # lower bounds for global effects
      
      emp_var <- var(svcov.object@z)
      
      max_dist <- max(as.matrix(dist(svcov.object@locs))) * 0.7
      
      limits <- rbind(mean.limits,
                      std.dev.max.effects,
                      scale.max.effects,
                      aniso.max.effects,
                      tilt.max.effects,
                      smooth.max.effects,
                      nugget.max.effects)
      
      tmp_par.pos <- getDesignMatrix(model.list = svcov.object@model.list,
                                     data = svcov.object@data)$par.pos
      
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
        
        names(nugget.max.effects) <- 'nugget.max.effects'
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
      
      return(list('theta_init' = theta_init,
                  'theta_upper' = theta_upper,
                  'theta_lower' = theta_lower))
      
    }
    
  } else{stop('svcov.object required')}
}


#' getHessian
#' @description returns the approximate (observed) Hesian (inverse of Fisher Information Matrix)
#' @usage getHessian(svcov.object, ncores = 2, 
#' eps = .Machine$double.eps^(1/4))
#' @param svcov.object a fitted svcov object
#' @param ncores number of cores used for the computation
#' @param eps ...
#' Fisher Information Matrix
#' @returns a symmetric matrix pxp of the approximated (observed) Hessian
#' @author Federico Blasi
getHessian <- function(svcov.object, ncores = parallel::detectCores() - 1, 
                       eps = .Machine$double.eps^(1/4)){
  
  if(svcov.object@type == 'dense'){
    
    f00 <- svcov.object@output$value
    
    p <- base::length(svcov.object@output$par)
    H <- base::matrix(NA, p, p)
    
    index_positions <- base::matrix(ncol = 2)
    
    for(jj in 1:(p)){
      for(ii in (jj):p){
        index_positions <- base::rbind(index_positions, c(jj, ii))
      }
    }
    
    index_positions <- index_positions[-1, ]
    
    Design_Matrix <- svcov::getDesignMatrix(model.list = svcov.object@model.list, 
                                            data = svcov.object@data)
    par.pos <- Design_Matrix$par.pos
    x_covariates <- Design_Matrix$model.matrix
    
    svcov_x_std <- svcov::getScale(x_covariates, 
                                   mean.vector = svcov.object@info$mean.vector,
                                   sd.vector = svcov.object@info$sd.vector)
    
    z <- svcov.object@z
    
    n <- dim(svcov.object@z)[1]
    
    x_covariates <- svcov_x_std$std.covs
    
    locs <- svcov.object@locs
    
    lambda <- svcov.object@info$lambda
    
    svcov.info <- svcov.object@info$smooth_limits
    
    pars <- svcov.object@output$par
    
    cl <- parallel::makeCluster(ncores)
    parallel::setDefaultCluster(cl = cl)
    parallel::clusterEvalQ(cl, library('svcov'))
    
    parallel::clusterExport(cl = cl, list("z", "x_covariates", "svcov.info",
                                          "locs", "n", "par.pos", "eps",
                                          "pars", "f00","lambda"),
                            envir = environment())
    
    vector_responses <- parallel::parApply(cl = cl, 
                                           index_positions, 
                                           MARGIN = 1, 
                                           FUN = function(x){
                                             
                                             t01 <- pars
                                             t10 <- pars
                                             t11 <- pars
                                             
                                             t01[x[1]] <- t01[x[1]] + eps
                                             t10[x[2]] <- t10[x[2]] + eps
                                             t11[x[1]] <- t11[x[1]] + eps
                                             t11[x[2]] <- t11[x[2]] + eps
                                             
                                             f01 <- svcov::GetNeg2loglikelihood(t01, par.pos = par.pos,
                                                                                    locs = locs,
                                                                                    x_covariates = x_covariates, 
                                                                                    smooth.limits = svcov.info, 
                                                                                    n = n, 
                                                                                    z = z,
                                                                                    lambda = lambda)
                                             
                                             f10 <- svcov::GetNeg2loglikelihood(t10, par.pos = par.pos,
                                                                                    locs = locs,
                                                                                    x_covariates = x_covariates, 
                                                                                    smooth.limits = svcov.info, 
                                                                                    n = n, 
                                                                                    z = z,
                                                                                    lambda = lambda)
                                             
                                             f11 <- svcov::GetNeg2loglikelihood(t11, par.pos = par.pos,
                                                                                    locs = locs,
                                                                                    x_covariates = x_covariates, 
                                                                                    smooth.limits = svcov.info, 
                                                                                    n = n, 
                                                                                    z = z, 
                                                                                    lambda = lambda)
                                             
                                             0.5 * ((f11 - f01 - f10 + f00) / (eps*eps))
                                           }); parallel::stopCluster(cl)
    
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
  
  if(svcov.object@type == 'sparse'){
    
    f00 <- svcov.object@output$value
    
    p <- base::length(svcov.object@output$par)
    H <- base::matrix(NA, p, p)
    
    index_positions <- base::matrix(ncol = 2)
    
    for(jj in 1:(p)){
      for(ii in (jj):p){
        index_positions <- base::rbind(index_positions, c(jj, ii))
      }
    }
    
    index_positions <- index_positions[-1, ]
    
    Design_Matrix <- svcov::getDesignMatrix(model.list = svcov.object@model.list, 
                                            data = svcov.object@data)
    par.pos <- Design_Matrix$par.pos
    x_covariates <- Design_Matrix$model.matrix
    
    svcov_x_std <- svcov::getScale(x_covariates, 
                                   mean.vector = svcov.object@info$mean.vector,
                                   sd.vector = svcov.object@info$sd.vector)
    
    z <- svcov.object@z
    n <- dim(svcov.object@z)[1]
    x_covariates <- svcov_x_std$std.covs
    
    locs <- svcov.object@locs
    
    lambda <- svcov.object@info$lambda
    
    svcov.info <- svcov.object@info$smooth_limits
    
    pars <- svcov.object@output$par
    
    ref_taper <- svcov.object@info$taper(
      spam::nearest.dist(svcov.object@locs, delta = svcov.object@info$delta, upper = NULL),
      theta = c(svcov.object@info$delta, 1)
    )
    
    cholS <- spam::chol.spam(ref_taper)

    cl <- parallel::makeCluster(ncores)
    parallel::setDefaultCluster(cl = cl)
    parallel::clusterEvalQ(cl, library('svcov'))
    
    parallel::clusterExport(cl = cl, list("z", "x_covariates", "svcov.info",
                                          "locs", "n", "par.pos", "eps",
                                          "pars", "f00","lambda","ref_taper","cholS"),
                            envir = environment())
    
    vector_responses <- parallel::parApply(cl = cl, 
                                           index_positions, 
                                           MARGIN = 1, 
                                           FUN = function(x){
                                             
                                             t01 <- pars
                                             t10 <- pars
                                             t11 <- pars
                                             
                                             t01[x[1]] <- t01[x[1]] + eps
                                             t10[x[2]] <- t10[x[2]] + eps
                                             t11[x[1]] <- t11[x[1]] + eps
                                             t11[x[2]] <- t11[x[2]] + eps
                                             
                                             f01 <- svcov::GetNeg2loglikelihoodTaper(theta = t01,
                                                                                          par.pos = par.pos,
                                                                                          ref_taper = ref_taper,
                                                                                          locs = locs,
                                                                                          x_covariates = x_covariates,
                                                                                          smooth.limits = svcov.info,
                                                                                          cholS = cholS,
                                                                                          z = z,
                                                                                          n = n,
                                                                                          lambda = lambda)
                                             
                                             f10 <- svcov::GetNeg2loglikelihoodTaper(theta = t10,
                                                                                          par.pos = par.pos,
                                                                                          ref_taper = ref_taper,
                                                                                          locs = locs,
                                                                                          x_covariates = x_covariates,
                                                                                          smooth.limits = svcov.info,
                                                                                          cholS = cholS,
                                                                                          z = z,
                                                                                          n = n,
                                                                                          lambda = lambda)
                                             
                                             f11 <- svcov::GetNeg2loglikelihoodTaper(theta = t11,
                                                                                          par.pos = par.pos,
                                                                                          ref_taper = ref_taper,
                                                                                          locs = locs,
                                                                                          x_covariates = x_covariates,
                                                                                          smooth.limits = svcov.info,
                                                                                          cholS = cholS,
                                                                                          z = z,
                                                                                          n = n,
                                                                                          lambda = lambda)
                                             
                                             0.5 * ((f11 - f01 - f10 + f00) / (eps*eps))
                                             
                                           }); parallel::stopCluster(cl)
    
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
