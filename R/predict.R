
#' Prediction for object class coco
#' @description Prediction for a fitted coco object.
#' @usage cocoPredict(coco.object, newdataset, newlocs, type = 'mean')
#' @param coco.object a coco object
#' @param newdataset a data.frame including covariates present in model.list
#' @param newlocs a matrix with locations related to newdataset
#' @param type whether "mean" or "pred", which gives a point prediction as well as the prediction uncertainty
#' @returns a list with mean predictions and CRPS values, latter if specified.
#' @author Federico Blasi
#' 
cocoPredict <- function(coco.object, 
                        newdataset,
                        newlocs,
                        type = "mean") {
  
  if (!('coco' %in% class(coco.object))){
    stop("coco class required.")
  }
  
  if (length(coco.object@output) == 0) {
    stop("object has not yet been fitted.")
  }
  
  .coco.check.newdataset(newdataset)
  .coco.check.newlocs(newlocs)
  .coco.check.type_pred(type)
  .coco.check.object(coco.object)
  
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
      smooth_limits = coco.object@info$smooth_limits
    )
    
    cov_pred <- cocons::cov_rns_pred(
      theta = adjusted_eff_values[-1], locs = coco.object@locs,
      locs_pred = as.matrix(newlocs),
      x_covariates = X_std$std.covs,
      x_covariates_pred = X_pred_std$std.covs,
      smooth_limits = coco.object@info$smooth_limits
    )
    
    inv_cov <- solve(observed_cov, t(cov_pred))
    
    # trend
    trend_pred <- c(X_pred_std$std.covs %*% adjusted_eff_values$mean)
    trendObs <- c(X_std$std.covs %*% adjusted_eff_values$mean)
    
    coco.resid <- coco.object@z[,1] - trendObs # Future updates add slice to pick from the samples
    
    # mean part
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
    taper_two@entries <- taper_two@entries * cocons::cov_rns_taper_optimized_range(
      theta = adjusted_eff_values,
      locs = coco.object@locs,
      x_covariates = X_std$std.covs,
      colindices = taper_two@colindices,
      rowpointers = taper_two@rowpointers,
      smooth_limits = coco.object@info$smooth_limits
    )
    
    pred_locs <- spam::nearest.dist(x = as.matrix(newlocs), y = coco.object@locs, delta = coco.object@info$delta)
    
    pred_taper <- coco.object@info$taper(pred_locs, theta = c(coco.object@info$delta, 1))
    
    rm(pred_locs)
    
    # C(preds, locs)
    pred_taper@entries <- pred_taper@entries * cocons::cov_rns_taper_optimized_predict_range(
      theta = adjusted_eff_values,
      locs = coco.object@locs,
      locs_pred = as.matrix(newlocs),
      x_covariates = X_std$std.covs,
      x_covariates_pred = X_pred_std$std.covs,
      colindices = pred_taper@colindices,
      rowpointers = pred_taper@rowpointers,
      smooth_limits = coco.object@info$smooth_limits
    )
    
    inv_cov <- spam::solve(taper_two, spam::t(pred_taper)) # consumes a lot of memory... more efficient would be solve(taper_two, coco.object@z) and then multiply... but for CRPS we need what it's implemented
    
    # trend
    trend_pred <- c(X_pred_std$std.covs %*% adjusted_eff_values$mean) # crossprod ?
    trend_obs <- c(X_std$std.covs %*% adjusted_eff_values$mean) # crossprod?
    coco.resid <- coco.object@z[,1,drop=FALSE] - trend_obs
    
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
