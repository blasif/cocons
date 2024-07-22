
#' Marginal simulation of Gaussian processes with nonstationary covariance function
#' @description simulates a Gaussian process with nonstationary covariance function
#' from an coco object.
#' 
#' @usage cocoSim(coco.object, pars, n, seed, standardize, 
#' type = 'classic', sim.type = NULL, cond.info = NULL)
#' @param coco.object a vector of length p, where p is the number of parameters for
#' @param pars a vector of length p, where p is the number of parameters for
#' each of the models
#' @param n number of realizations to simulate
#' @param seed a list detailing in which position of each aspect the elements
#' of theta should be placed. Expected to be output of getDesignMatrix
#' @param seed seed number. defalt set to NULL.
#' @param standardize logical argument describing whether provided covariates
#' should be standardize (TRUE) or not (FALSE). By default set to TRUE
#' @param type wether parameters are related to a classical parameterization ('classic') or
#' a difference parameterization 'diff' . Default set to 'classic'.
#' @param sim.type if set 'cond' then a conditional simulation takes place.
#' @param cond.info a list containing information to perform a conditional simulation.
#' @returns a list() with the structure needed to fit cov.rns functions
#' @author Federico Blasi
#' 
cocoSim <- function(coco.object, 
                     pars,
                     n = 1,
                     seed = NULL, 
                     standardize = TRUE, 
                     type = 'classic',
                     sim.type = NULL,
                     cond.info = NULL){
  
  # add a check to test whether length of pars match model specification
  
  if(!(coco.object@type %in% c('dense', 'sparse'))){stop('')} # should not be necessary if using a 
  # valid coco object
  
  if(coco.object@type == 'dense'){
    
    if(!is.null(sim.type)){
      if(sim.type == 'cond'){
        
        # check object fitted
        
        # test new
        
        std_coco <- getScale(coco.object)$std.covs
        
        test_here <- getDesignMatrix(coco.object@model.list, data = cond.info$newdataset)
        
        std_pred <- getScale(test_here$model.matrix, mean.vector = coco.object@info$mean.vector,
                             sd.vector = coco.object@info$sd.vector)$std.covs
        
        to_pass <- cocons::getModelLists(coco.object@output$par, 
                                        par.pos = test_here$par.pos, 
                                        type = 'diff')
        
        covmat <- cocons::cov_rns(theta = to_pass[-1], 
                                 locs = coco.object@locs,
                                 x_covariates = std_coco,
                                 smooth_limits = coco.object@info$smooth_limits)
        
        covmat_pred <- cocons::cov_rns_pred(theta = to_pass[-1], 
                                           locs = coco.object@locs,
                                           locs_pred = as.matrix(cond.info$newlocs),
                                           x_covariates = std_coco,
                                           x_covariates_pred = std_pred,
                                           smooth_limits = coco.object@info$smooth_limits)
        
        covmat_unobs <- cocons::cov_rns(theta = to_pass[-1], 
                                       locs = as.matrix(cond.info$newdataset),
                                       x_covariates = std_pred,
                                       smooth_limits = coco.object@info$smooth_limits)
        
        part_b <- covmat_pred %*% solve(covmat) %*% t(covmat_pred)
        
        conditional_covariance <-  covmat_unobs - part_b
        
        L <- chol(conditional_covariance)
        
        iiderrors <- replicate(n, expr = stats::rnorm(dim(cond.info$newlocs)[1], mean = 0, sd = 1))
        
        step_one <- cocons::cocoPredict(coco.object, 
                                        newdataset = cond.info$newdataset, 
                                        newlocs = as.matrix(cond.info$newlocs), 
                                        type = 'mean')
        
        tmp_mu <- step_one$trend + step_one$mean
        
        return(sweep(t(iiderrors) %*% L, 2, tmp_mu, "+"))
        
        if(F){
          # end test new
          
          step_one <- cocons::cocoPredict(coco.object, 
                                          newdataset = cond.info$newdataset, 
                                          newlocs = as.matrix(cond.info$newlocs), 
                                          type = 'mean')
          
          step_two <- cocons::cocoSim(coco.object, pars = pars, n = n, seed = seed,
                                      standardize = standardize, type = type)
          
          tmp_coco_object <- coco(type = 'dense',
                                    locs = as.matrix(cond.info$newlocs),
                                    data = cond.info$newdataset,
                                    z = numeric(length = dim(cond.info$newdataset)[1]),
                                    model.list = coco.object@model.list,
                                    info = coco.object@info)
          
          step_three <- cocons::cocoSim(tmp_coco_object, pars = pars, n = n, seed = seed, 
                                        standardize = standardize, type = type)
          
          matrix_return <- matrix(NA, nrow = dim(step_three)[1], ncol = dim(step_three)[2])
          
          coco.object_four <- coco.object
          
          for(ii in 1:n){
            
            coco.object_four@z <- c(step_two[ii, ])
            
            step_four <- cocons::cocoPredict(coco.object_four, 
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
      
      coco_items <- getDesignMatrix(model.list = coco.object@model.list, 
                                     data = coco.object@data)
      
      if(standardize){
        std_coco <- cocons::getScale(coco_items$model.matrix)
      } else{
        std_coco <- cocons::getScale(coco_items$model.matrix,
                                     mean.vector = rep(0, dim(coco_items$model.matrix)[2]),
                                     sd.vector = rep(1, dim(coco_items$model.matrix)[2]))
      }
      
      theta_to_fit <- cocons::getModelLists(pars,
                                           par.pos = coco_items$par.pos, 
                                           type = type)
      
      if(type == 'classic'){
        
        if(!is.formula(coco.object@model.list$smooth)){
          
          theta_to_fit$smooth[1] <- log(coco.object@info$smooth_limits[1])
          
        }
        
        covmat <- cocons::cov_rns_classic(theta = theta_to_fit[-1], 
                                         locs = coco.object@locs,
                                         x_covariates = std_coco$std.covs) 
        
      }
      
      if(type == 'diff'){
        covmat <- cocons::cov_rns(theta = theta_to_fit[-1], 
                                 locs = coco.object@locs,
                                 x_covariates = std_coco$std.covs,
                                 smooth_limits = coco.object@info$smooth_limits) 
      }
      
      cholS <- base::chol(covmat)
      
      if(exists("seed")){
        set.seed(seed)
      }
      
      iiderrors <- replicate(n, expr = stats::rnorm(dim(coco.object@data)[1], mean = 0, sd = 1))
      
      tmp_mu <- std_coco$std.covs %*% theta_to_fit$mean
      
      return(sweep(t(iiderrors) %*% cholS, 2, tmp_mu, "+"))
      
    }
  } 
  
  if(coco.object@type == 'sparse'){
    
    coco_items <- getDesignMatrix(model.list = coco.object@model.list, data = coco.object@data)
    
    if(standardize){
      std_coco <- cocons::getScale(coco_items$model.matrix)
    } else{
      std_coco <- cocons::getScale(coco_items$model.matrix,
                                   mean.vector = rep(0, dim(coco_items$model.matrix)[2]),
                                   sd.vector = rep(1, dim(coco_items$model.matrix)[2]))
    }
    
    theta_to_fit <- cocons::getModelLists(pars,
                                         par.pos = coco_items$par.pos, 
                                         type = type)
    
    ref_taper <- coco.object@info$taper(
      spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL), 
      theta = c(coco.object@info$delta, 1)
    )
    
    ref_taper@entries <- ref_taper@entries * cov_rns_taper_optimized_range(theta = theta_to_fit[-1], 
                                                                           locs = coco.object@locs, 
                                                                           x_covariates =  std_coco$std.covs, 
                                                                           colindices = ref_taper@colindices, 
                                                                           rowpointers = ref_taper@rowpointers,
                                                                           smooth_limits =  coco.object@info$smooth_limits)
    cholS <- spam::chol(ref_taper) # check rmvnorm.spam
    iord <- spam::ordering(cholS, inv = TRUE)
    cholS <- spam::as.spam(cholS)
    
    if(exists("seed")){
      set.seed(seed)
    }
    
    iiderrors <- replicate(n, expr = stats::rnorm(dim(coco.object@data)[1], mean = 0, sd = 1))
    
    tmp_mu <- std_coco$std.covs %*% theta_to_fit$mean
    
    return(sweep(as.matrix(t(iiderrors) %*% cholS)[,iord, drop = F], 2, tmp_mu, "+"))
    
  }
  
}
