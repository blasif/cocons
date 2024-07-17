
#' Marginal simulation of Gaussian Processes with nonstationary covariance function
#' @description simulates a Gaussian Process with nonstationary covariance function
#' from an svcov object.
#' 
#' @usage svcovSim(svcov.object, pars, n, seed, standardize, type = 'classic')
#' @param svcov.object a vector of length p, where p is the number of parameters for
#' @param pars a vector of length p, where p is the number of parameters for
#' each of the models
#' @param n number of realizations to simulate
#' @param seed a list detailing in which position of each aspect the elements
#' of theta should be placed. Expected to be output of getDesignMatrix
#' @param cond.info a list containing information to perform a conditional simulation.
#' @param standardize logical argument describing whether provided covariates
#' should be standardize (TRUE) or not (FALSE). By default set to TRUE
#' @param seed seed number. defalt set to NULL.
#' @param type wether parameters are related to a classical parameterization ('classic') or
#' a difference parameterization 'diff' . Default set to 'classic'.
#' @param sim.type if set 'cond' then a conditional simulation takes place.
#' @returns a list() with the structure needed to fit cov.rns functions
#' @examples 
#' model.list <- list("std.dev" = as.formula(" ~ 1 + lati_s * long_s"),
#'                       "scale" = as.formula(" ~ 1 + elev_s"),
#'                       "aniso" = as.formula(" ~ 1 + elev_s"),
#'                       "tilt" = as.formula(" ~ 1 + elev_s"),
#'                       "smooth" = as.formula(" ~ 1"),
#'                       "nugget" = -Inf)
#' @author Federico Blasi
#' 
svcovSim <- function(svcov.object, 
                     pars,
                     n = 1,
                     seed = NULL, 
                     standardize = TRUE, 
                     type = 'classic',
                     sim.type = NULL,
                     cond.info = NULL){
  
  # add a check to test whether length of pars match model specification
  
  if(!(svcov.object@type %in% c('dense', 'sparse'))){stop('')} # should not be necessary if using a 
  # valid svcov object
  
  if(svcov.object@type == 'dense'){
    
    if(!is.null(sim.type)){
      if(sim.type == 'cond'){
        
        # check object fitted
        
        # test new
        
        std_svcov <- getScale(svcov.object)$std.covs
        
        test_here <- getDesignMatrix(svcov.object@model.list, data = cond.info$newdataset)
        
        std_pred <- getScale(test_here$model.matrix, mean.vector = svcov.object@info$mean.vector,
                             sd.vector = svcov.object@info$sd.vector)$std.covs
        
        to_pass <- svcov::getModelLists(svcov.object@output$par, 
                                        par.pos = test_here$par.pos, 
                                        type = 'diff')
        
        covmat <- svcov::cov_rns(theta = to_pass[-1], 
                                 locs = svcov.object@locs,
                                 x_covariates = std_svcov,
                                 smooth_limits = svcov.object@info$smooth_limits)
        
        covmat_pred <- svcov::cov_rns_pred(theta = to_pass[-1], 
                                           locs = svcov.object@locs,
                                           locs_pred = as.matrix(cond.info$newlocs),
                                           x_covariates = std_svcov,
                                           x_covariates_pred = std_pred,
                                           smooth_limits = svcov.object@info$smooth_limits)
        
        covmat_unobs <- svcov::cov_rns(theta = to_pass[-1], 
                                       locs = as.matrix(cond.info$newdataset),
                                       x_covariates = std_pred,
                                       smooth_limits = svcov.object@info$smooth_limits)
        
        part_b <- covmat_pred %*% solve(covmat) %*% t(covmat_pred)
        
        conditional_covariance <-  covmat_unobs - part_b
        
        L <- chol(conditional_covariance)
        
        iiderrors <- replicate(n, expr = stats::rnorm(dim(cond.info$newlocs)[1], mean = 0, sd = 1))
        
        step_one <- svcov::svcovPredict(svcov.object, 
                                        newdataset = cond.info$newdataset, 
                                        newlocs = as.matrix(cond.info$newlocs), 
                                        type = 'mean')
        
        tmp_mu <- step_one$trend + step_one$mean
        
        return(sweep(t(iiderrors) %*% L, 2, tmp_mu, "+"))
        
        if(F){
          # end test new
          
          step_one <- svcov::svcovPredict(svcov.object, 
                                          newdataset = cond.info$newdataset, 
                                          newlocs = as.matrix(cond.info$newlocs), 
                                          type = 'mean')
          
          step_two <- svcov::svcovSim(svcov.object, pars = pars, n = n, seed = seed,
                                      standardize = standardize, type = type)
          
          tmp_svcov_object <- svcov(type = 'dense',
                                    locs = as.matrix(cond.info$newlocs),
                                    data = cond.info$newdataset,
                                    z = numeric(length = dim(cond.info$newdataset)[1]),
                                    model.list = svcov.object@model.list,
                                    info = svcov.object@info)
          
          step_three <- svcov::svcovSim(tmp_svcov_object, pars = pars, n = n, seed = seed, 
                                        standardize = standardize, type = type)
          
          matrix_return <- matrix(NA, nrow = dim(step_three)[1], ncol = dim(step_three)[2])
          
          svcov.object_four <- svcov.object
          
          for(ii in 1:n){
            
            svcov.object_four@z <- c(step_two[ii, ])
            
            step_four <- svcov::svcovPredict(svcov.object_four, 
                                             newdataset = cond.info$newdataset, 
                                             newlocs = as.matrix(cond.info$newlocs),
                                             type = 'mean')
            
            matrix_return[ii, ] <- c(step_one$mean + step_one$trend) + c(step_three[ii, ]) - 
              c(step_four$trend + step_four$mean) # trend should not be part of the simulation here I think!
            
          }
          
          return(matrix_return)
        }
      } 
    } else{
      
      svcov_items <- getDesignMatrix(model.list = svcov.object@model.list, 
                                     data = svcov.object@data)
      
      if(standardize){
        std_svcov <- svcov::getScale(svcov_items$model.matrix)
      } else{
        std_svcov <- svcov::getScale(svcov_items$model.matrix,
                                     mean.vector = rep(0, dim(svcov_items$model.matrix)[2]),
                                     sd.vector = rep(1, dim(svcov_items$model.matrix)[2]))
      }
      
      theta_to_fit <- svcov::getModelLists(pars,
                                           par.pos = svcov_items$par.pos, 
                                           type = type)
      
      if(type == 'classic'){
        
        if(!is.formula(svcov.object@model.list$smooth)){
          
          theta_to_fit$smooth[1] <- log(svcov.object@info$smooth_limits[1])
          
        }
        
        covmat <- svcov::cov_rns_classic(theta = theta_to_fit[-1], 
                                         locs = svcov.object@locs,
                                         x_covariates = std_svcov$std.covs) 
        
      }
      
      if(type == 'diff'){
        covmat <- svcov::cov_rns(theta = theta_to_fit[-1], 
                                 locs = svcov.object@locs,
                                 x_covariates = std_svcov$std.covs,
                                 smooth_limits = svcov.object@info$smooth_limits) 
      }
      
      cholS <- base::chol(covmat)
      
      if(exists("seed")){
        set.seed(seed)
      }
      
      iiderrors <- replicate(n, expr = stats::rnorm(dim(svcov.object@data)[1], mean = 0, sd = 1))
      
      tmp_mu <- std_svcov$std.covs %*% theta_to_fit$mean
      
      return(sweep(t(iiderrors) %*% cholS, 2, tmp_mu, "+"))
      
    }
  } 
  
  if(svcov.object@type == 'sparse'){
    
    svcov_items <- getDesignMatrix(model.list = svcov.object@model.list, data = svcov.object@data)
    
    if(standardize){
      std_svcov <- svcov::getScale(svcov_items$model.matrix)
    } else{
      std_svcov <- svcov::getScale(svcov_items$model.matrix,
                                   mean.vector = rep(0, dim(svcov_items$model.matrix)[2]),
                                   sd.vector = rep(1, dim(svcov_items$model.matrix)[2]))
    }
    
    theta_to_fit <- svcov::getModelLists(pars,
                                         par.pos = svcov_items$par.pos, 
                                         type = type)
    
    ref_taper <- svcov.object@info$taper(
      spam::nearest.dist(svcov.object@locs, delta = svcov.object@info$delta, upper = NULL), 
      theta = c(svcov.object@info$delta, 1)
    )
    
    ref_taper@entries <- ref_taper@entries * cov_rns_taper_optimized_range(theta = theta_to_fit[-1], 
                                                                           locs = svcov.object@locs, 
                                                                           x_covariates =  std_svcov$std.covs, 
                                                                           colindices = ref_taper@colindices, 
                                                                           rowpointers = ref_taper@rowpointers,
                                                                           smooth_limits =  svcov.object@info$smooth_limits)
    cholS <- spam::chol(ref_taper) # check rmvnorm.spam
    iord <- spam::ordering(cholS, inv = TRUE)
    cholS <- spam::as.spam(cholS)
    
    if(exists("seed")){
      set.seed(seed)
    }
    
    iiderrors <- replicate(n, expr = stats::rnorm(dim(svcov.object@data)[1], mean = 0, sd = 1))
    
    tmp_mu <- std_svcov$std.covs %*% theta_to_fit$mean
    
    return(sweep(as.matrix(t(iiderrors) %*% cholS)[,iord, drop = F], 2, tmp_mu, "+"))
    
  }
  
}
