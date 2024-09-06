
#' Prediction routines for nonstationary spatial models
#' @description Computes the point predictions and standard errors based on conditional Gaussian distributions.
#' @usage cocoPredict(coco.object, newdataset, newlocs, type = 'mean', ...)
#' @param coco.object (\code{S4}) a fitted \link{coco} object.
#' @param newdataset (\code{data.frame}) a data.frame containing covariates present in model.list at prediction locations.
#' @param newlocs (\code{matrix}) a matrix with locations related to prediction locations, matching indexing of \code{newdataset}.
#' @param type (\code{character}) whether \code{"mean"} or \code{"pred"}, which gives a point prediction for the former, 
#' as well as of point prediction and standard errors for the latter.
#' @param ... (\code{character}) when coco.object has multiple realizations, specifying \code{"index.pred"} specifying which column of 
#' \code{coco.object@z} should be used to perform predictions.
#' @returns (\code{list}) a list with the conditional mean, splitted into the systematic large-scale variability \code{trend}, 
#' and due to stochastic \code{mean}, as well as standard errors \code{"sd.pred"} if \code{"pred"} is specified.
#' @author Federico Blasi
#' @examples
#' \dontrun{
#'  
#' model.list <- list('mean' = 0,
#' 'std.dev' = formula( ~ 1 + cov_x + cov_y),
#' 'scale' = formula( ~ 1 + cov_x + cov_y),
#' 'aniso' = 0,
#' 'tilt' = 0,
#' 'smooth' = 3/2,
#' 'nugget' = -Inf)
#' 
#' coco_object <- coco(type = 'dense',
#' data = holes[[1]][1:100, ],
#' locs = as.matrix(holes[[1]][1:100, 1:2]),
#' z = holes[[1]][1:100, ]$z,
#' model.list = model.list)
#' 
#' optim_coco <- cocoOptim(coco_object,
#' boundaries = getBoundaries(coco_object,
#' lower.value = -3, 3))
#' 
#' coco_preds <- cocoPredict(optim_coco, newdataset = holes[[2]],
#' newlocs = as.matrix(holes[[2]][, 1:2]),
#' type = "pred")
#' 
#' coco_preds
#' 
#' par(mfrow = c(1, 2))
#' 
#' fields::quilt.plot(main = "mean", holes[[2]][, 1:2], 
#' coco_preds$mean, xlim = c(-1, 1), ylim = c(-1, 1))
#' fields::quilt.plot(main = "se", holes[[2]][, 1:2], 
#' coco_preds$sd.pred, xlim = c(-1, 1), ylim = c(-1, 1))
#' 
#' # Re-do it without considering cov_x and cov_y in the std.dev and scale and compare. 
#' 
#' }
#' 
cocoPredict <- function(coco.object, 
                        newdataset,
                        newlocs,
                        type = "mean",
                        ...) {
  
  .cocons.check.coco(coco.object)

  if (length(coco.object@output) == 0) {
    stop("coco object has not yet been fitted.")
  }
  
  if(dim(coco.object@z)[2] == 1){
    index.pred <- 1
  } else{
    
    if(!exists("index.pred")){
      index.pred <- 1
    } else{
      if(index.pred > dim(coco.object@z)[2]){
        stop("index.pred is larger than dim(coco.object@z)[2]")
      }
      }
  }
  
  .cocons.check.newdataset(newdataset, coco.object)
  .cocons.check.newlocs(newlocs)
  .cocons.check.type_pred(type)

  # add check on the names of newdataset names and model.list
  # add check type
  
  if (coco.object@type == "dense") {
    
    tmp_matrix <- cocons::getDesignMatrix(model.list = coco.object@model.list, data = coco.object@data)
    
    adjusted_eff_values <- cocons::getModelLists(coco.object@output$par, 
                                                par.pos = tmp_matrix$par.pos, 
                                                type = "diff")
    
    X_std <- cocons::getScale(tmp_matrix$model.matrix,
                             mean.vector = coco.object@info$mean.vector,
                             sd.vector = coco.object@info$sd.vector
    )
    
    tmp_matrix_pred <- cocons::getDesignMatrix(
      model.list = coco.object@model.list,
      data = newdataset
    )
    
    X_pred_std <- cocons::getScale(tmp_matrix_pred$model.matrix,
                                  mean.vector = coco.object@info$mean.vector,
                                  sd.vector = coco.object@info$sd.vector
    )
    
    observed_cov <- cocons::cov_rns(
      theta = adjusted_eff_values[-1], locs = coco.object@locs,
      x_covariates = X_std$std.covs,
      smooth_limits = coco.object@info$smooth.limits
    )
    
    cov_pred <- cocons::cov_rns_pred(
      theta = adjusted_eff_values[-1], locs = coco.object@locs,
      locs_pred = as.matrix(newlocs),
      x_covariates = X_std$std.covs,
      x_covariates_pred = X_pred_std$std.covs,
      smooth_limits = coco.object@info$smooth.limits
    )
    
    inv_cov <- solve(observed_cov, t(cov_pred))
    
    # trend
    trend_pred <- c(X_pred_std$std.covs %*% adjusted_eff_values$mean)
    trendObs <- c(X_std$std.covs %*% adjusted_eff_values$mean)
    
    coco.resid <- coco.object@z[,index.pred] - trendObs
    
    # spatial mean
    mean_part <- c(crossprod(coco.resid, inv_cov))
    
    if (type == "mean") {
      return(list(
        "trend" = trend_pred,
        "mean" = mean_part
      ))
    }
    
    uncertainty_some <- 1 / exp(-X_pred_std$std.covs %*% adjusted_eff_values$std.dev) +
      exp(X_pred_std$std.covs %*% adjusted_eff_values$nugget)

    if (type == "pred") {
      
      vector_tmp_z <- numeric(dim(newlocs)[1])
      
      for (ii in 1:dim(newlocs)[1]) {
        uncertainty_some[ii] <- uncertainty_some[ii] - cov_pred[ii, , drop = FALSE] %*% inv_cov[, ii, drop = FALSE]
        if(abs(uncertainty_some[ii]) < 1e-10) uncertainty_some[ii] <- abs(uncertainty_some[ii]) # rounding errors
      }
      
      return(
        list(
          "trend" = trend_pred,
          "mean" = mean_part,
          "sd.pred" = c(sqrt(uncertainty_some))
        )
      )
      
    }
  }
  
  if (coco.object@type == "sparse") {
    
    tmp_matrix <- cocons::getDesignMatrix(
      model.list = coco.object@model.list,
      data = coco.object@data
    )
    
    adjusted_eff_values <- cocons::getModelLists(
      theta = coco.object@output$par,
      par.pos = tmp_matrix$par.pos, type = "diff"
    )
    
    tmp_matrix_pred <- cocons::getDesignMatrix(
      model.list = coco.object@model.list,
      data = newdataset
    )
    
    X_std <- cocons::getScale(tmp_matrix$model.matrix,
                             mean.vector = coco.object@info$mean.vector,
                             sd.vector = coco.object@info$sd.vector
    )
    
    X_pred_std <- cocons::getScale(tmp_matrix_pred$model.matrix,
                                  mean.vector = coco.object@info$mean.vector,
                                  sd.vector = coco.object@info$sd.vector
    )
    
    ###
    
    distmat <- spam::nearest.dist(coco.object@locs, delta = coco.object@info$delta, upper = NULL)
    
    taper_two <- coco.object@info$taper(distmat, theta = c(coco.object@info$delta, 1))
    
    # C(locs,locs)
    taper_two@entries <- taper_two@entries * cocons::cov_rns_taper(
      theta = adjusted_eff_values,
      locs = coco.object@locs,
      x_covariates = X_std$std.covs,
      colindices = taper_two@colindices,
      rowpointers = taper_two@rowpointers,
      smooth_limits = coco.object@info$smooth.limits
    )
    
    pred_locs <- spam::nearest.dist(x = as.matrix(newlocs), y = coco.object@locs, delta = coco.object@info$delta)
    
    pred_taper <- coco.object@info$taper(pred_locs, theta = c(coco.object@info$delta, 1))
    
    rm(pred_locs)
    
    # C(preds, locs)
    pred_taper@entries <- pred_taper@entries * cocons::cov_rns_taper_pred(
      theta = adjusted_eff_values,
      locs = coco.object@locs,
      locs_pred = as.matrix(newlocs),
      x_covariates = X_std$std.covs,
      x_covariates_pred = X_pred_std$std.covs,
      colindices = pred_taper@colindices,
      rowpointers = pred_taper@rowpointers,
      smooth_limits = coco.object@info$smooth.limits
    )
    
    inv_cov <- spam::solve(taper_two, spam::t(pred_taper)) # memory intensive
    
    # trend
    trend_pred <- c(X_pred_std$std.covs %*% adjusted_eff_values$mean) # crossprod ?
    trend_obs <- c(X_std$std.covs %*% adjusted_eff_values$mean) # crossprod?
    coco.resid <- coco.object@z[,index.pred] - trend_obs
    
    # mean part
    
    mean_part <- c(crossprod(coco.resid, inv_cov))
    
    if (type == "mean") {
      return(list(
        "trend" = trend_pred,
        "mean" = mean_part
      ))
    }
    
    uncertainty_some <- 1 / exp(-X_pred_std$std.covs %*% adjusted_eff_values$std.dev) + exp(X_pred_std$std.covs %*% adjusted_eff_values$nugget)
    
    if (type == "pred") {
      
      for (ii in 1:dim(newlocs)[1]) {
        uncertainty_some[ii] <- uncertainty_some[ii] - pred_taper[ii, , drop = FALSE] %*% inv_cov[, ii, drop = FALSE]
        if(abs(uncertainty_some[ii]) < 1e-10) uncertainty_some[ii] <- abs(uncertainty_some[ii]) # rounding errors
      }
      
      return(list(
        "trend" = trend_pred,
        "mean" = mean_part,
        "sd.pred" = c(sqrt(uncertainty_some))
        )
      )
      }
    }
}
