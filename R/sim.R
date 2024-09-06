
#' Marginal and conditional simulation of nonstationary Gaussian process
#' @description draw realizations of nonstationary Gaussian processes with covariate-based covariance functions.
#' @details \code{'cond'} sim.type requires specifying in \code{'cond.info'} a list with \code{'newdataset'} a data.frame containing covariates present in model.list at simulation locations, 
#' and \code{'newlocs'} a matrix with locations related to the simulation locations, matching indexing of \code{'newdataset'}. 
#' 
#' \code{type = 'classic'} assumes a simpler parameterization for the covariance function, assuming log-parameterizations for \code{'std.dev'}, \code{'scale'}, and \code{'smooth'}.
#' 
#' @usage cocoSim(coco.object, pars, n, seed, standardize, 
#' type = 'classic', sim.type = NULL, cond.info = NULL)
#' @param coco.object (\code{S4}) a \link{coco} object.
#' @param pars (\code{numeric vector}) a vector of parameters values related to \code{model.list}.
#' @param n (\code{integer}) number of realizations to simulate.
#' @param seed (\code{integer or NULL}) seed number. default set to NULL.
#' @param standardize (\code{TRUE/FALSE}) logical argument describing whether provided covariates 
#' should be standardize (TRUE) or not (FALSE). By default set to TRUE.
#' @param type (\code{character}) whether parameters are related to a classical parameterization ('classic') or
#' a difference parameterization \code{'diff'}. Default set to \code{'classic'}. For \code{'sparse'} coco objects, only \code{'diff'} is available.
#' @param sim.type (\code{character}) if set \code{'cond'} then a conditional simulation takes place.
#' @param cond.info (\code{list}) a list containing added information to perform conditional simulation.
#' @returns (\code{matrix}) a matrix n x dim(data)\[1\].
#' @author Federico Blasi
#' @seealso [coco()]
#' @examples
#' \dontrun{
#' 
#' model.list <- list('mean' = 0,
#'                    'std.dev' = formula( ~ 1 + cov_x + cov_y),
#'                    'scale' = formula( ~ 1 + cov_x + cov_y),
#'                    'aniso' = 0,
#'                    'tilt' = 0,
#'                    'smooth' = 0.5,
#'                    'nugget' = -Inf)
#'                    
#' coco_object <- coco(type = 'dense',
#'                     data = holes[[1]][1:1000,],
#'                     locs = as.matrix(holes[[1]][1:1000,1:2]),
#'                     z = holes[[1]][1:1000,]$z,
#'                     model.list = model.list)
#'                     
#' coco_sim <- cocoSim(coco.object = coco_object,
#'             pars = c(0,0.25,0.25,  # pars related to std.dev
#'             log(0.25),1,-1),       # pars related to scale
#'             n = 1, 
#'             standardize = TRUE) 
#'
#' fields::quilt.plot(coco_object@locs,coco_sim)             
#' }
#' 
cocoSim <- function(coco.object, 
                     pars,
                     n = 1,
                     seed = NULL, 
                     standardize = TRUE, 
                     type = "classic",
                     sim.type = NULL,
                     cond.info = NULL){
  
  # add a check to test whether length of pars match model specification
  
  .cocons.check.type(coco.object@type)
  
  if(coco.object@type == "dense"){
    
    if(!is.null(sim.type)){
      if(sim.type == "cond"){
        
        # check object fitted
        
        std_coco <- getScale(coco.object)$std.covs
        
        test_here <- getDesignMatrix(coco.object@model.list, data = cond.info$newdataset)
        
        std_pred <- getScale(test_here$model.matrix, mean.vector = coco.object@info$mean.vector,
                             sd.vector = coco.object@info$sd.vector)$std.covs
        
        to_pass <- cocons::getModelLists(coco.object@output$par, 
                                        par.pos = test_here$par.pos, 
                                        type = "diff")
        
        covmat <- cocons::cov_rns(theta = to_pass[-1], 
                                 locs = coco.object@locs,
                                 x_covariates = std_coco,
                                 smooth_limits = coco.object@info$smooth.limits)
        
        covmat_pred <- cocons::cov_rns_pred(theta = to_pass[-1], 
                                           locs = coco.object@locs,
                                           locs_pred = as.matrix(cond.info$newlocs),
                                           x_covariates = std_coco,
                                           x_covariates_pred = std_pred,
                                           smooth_limits = coco.object@info$smooth.limits)
        
        covmat_unobs <- cocons::cov_rns(theta = to_pass[-1], 
                                       locs = as.matrix(cond.info$newdataset),
                                       x_covariates = std_pred,
                                       smooth_limits = coco.object@info$smooth.limits)
        
        part_b <- covmat_pred %*% solve(covmat) %*% t(covmat_pred)
        
        conditional_covariance <-  covmat_unobs - part_b
        
        L <- chol(conditional_covariance)
        
        iiderrors <- replicate(n, expr = stats::rnorm(dim(cond.info$newlocs)[1], mean = 0, sd = 1))
        
        step_one <- cocons::cocoPredict(coco.object, 
                                        newdataset = cond.info$newdataset, 
                                        newlocs = as.matrix(cond.info$newlocs), 
                                        type = "mean")
        
        tmp_mu <- step_one$trend + step_one$mean
        
        return(sweep(t(iiderrors) %*% L, 2, tmp_mu, "+"))
        
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
      
      if(type == "classic"){
        
        if(!is.formula(coco.object@model.list$smooth)){
          
          theta_to_fit$smooth[1] <- log(coco.object@info$smooth.limits[1])
          
        }
        
        covmat <- cocons::cov_rns_classic(theta = theta_to_fit[-1], 
                                         locs = coco.object@locs,
                                         x_covariates = std_coco$std.covs) 
        
      }
      
      if(type == "diff"){
        covmat <- cocons::cov_rns(theta = theta_to_fit[-1], 
                                 locs = coco.object@locs,
                                 x_covariates = std_coco$std.covs,
                                 smooth_limits = coco.object@info$smooth.limits) 
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
  
  if(coco.object@type == "sparse"){
    
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
    
    ref_taper@entries <- ref_taper@entries * cov_rns_taper(theta = theta_to_fit[-1], 
                                                                           locs = coco.object@locs, 
                                                                           x_covariates =  std_coco$std.covs, 
                                                                           colindices = ref_taper@colindices, 
                                                                           rowpointers = ref_taper@rowpointers,
                                                                           smooth_limits =  coco.object@info$smooth.limits)
    cholS <- spam::chol(ref_taper)
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
