
#' Marginal and conditional simulation of nonstationary Gaussian processes
#' @description draw realizations of stationary and nonstationary Gaussian processes with covariate-based covariance functions.
#' @details 
#' The argument \code{sim.type = 'cond'} specifies a conditional simulation, requiring \code{cond.info} to be provided. 
#' \code{cond.info} is a list including \code{newdataset}, a data.frame containing covariates present in \code{model.list} at the simulation locations, and \code{newlocs}, 
#' a matrix specifying the locations corresponding to the simulation, with indexing that matches \code{newdataset}.
#' 
#' The argument \code{type = 'classic'} assumes a simplified parameterization for the covariance function, with log-parameterizations applied to the parameters \code{std.dev}, 
#' \code{scale}, and \code{smooth}. 
#' 
#' @usage cocoSim(coco.object, pars, n, seed, standardize, 
#' type = 'classic', sim.type = NULL, cond.info = NULL)
#' @param coco.object (\code{S4}) A \link{coco} object.
#' @param pars (\code{numeric vector} or NULL) A vector of parameter values associated with \code{model.list}. 
#' If coco.object is a fitted object, and pars is \code{NULL}, it get pars from coco.object\@output$pars (and also sets 'type' to 'diff').
#' @param n (\code{integer}) Number of realizations to simulate.
#' @param seed (\code{integer or NULL}) Seed for random number generation. Defaults to NULL.
#' @param standardize (\code{logical}) Indicates whether the provided covariates should be standardized (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{TRUE}.
#' @param type (\code{character}) Specifies whether the parameters follow a classical parameterization (\code{'classic'}) or a difference parameterization (\code{'diff'}). Defaults to \code{'classic'}. For sparse \code{coco} objects, only \code{'diff'} is allowed.
#' @param sim.type (\code{character}) If set to \code{'cond'}, a conditional simulation is performed.
#' @param cond.info (\code{list}) A list containing additional information required for conditional simulation.
#' @returns (\code{matrix}) a matrix dim(data)\[1\] x n.
#' @author Federico Blasi
#' @seealso \link{coco}
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
                     pars = NULL,
                     n = 1,
                     seed = NULL, 
                     standardize = TRUE, 
                     type = "classic",
                     sim.type = NULL,
                     cond.info = NULL){
  
  .cocons.check.pars(coco.object,pars)
  .cocons.check.type(coco.object@type)
  
  if(is.null(pars)){
    pars <- coco.object@output$par
    type <- "diff"
  }
  
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
        
        
        
        L <- base::chol(covmat_unobs - covmat_pred %*% solve(covmat, t(covmat_pred))) 
        
        if(exists("seed")){
          set.seed(seed)
        }
        
        iiderrors <- replicate(n, expr = stats::rnorm(dim(cond.info$newlocs)[1], mean = 0, sd = 1))
        
        step_one <- cocons::cocoPredict(coco.object, 
                                        newdataset = cond.info$newdataset, 
                                        newlocs = as.matrix(cond.info$newlocs), 
                                        type = "mean")
        
        tmp_mu <- step_one$systematic + step_one$stochastic
        
        return(t(sweep(t(iiderrors) %*% L, 2, tmp_mu, "+")))
        
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
      
      if(!is.formula(coco.object@model.list$smooth)){
        
        theta_to_fit$smooth[1] <- log(coco.object@info$smooth.limits[1])
        
      }
      
      if(type == "classic"){
        
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
      
      return(t(sweep(t(iiderrors) %*% cholS, 2, tmp_mu, "+")))
      
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
    
    return(t(sweep(as.matrix(t(iiderrors) %*% cholS)[,iord, drop = F], 2, tmp_mu, "+")))
    
  }
  
}
