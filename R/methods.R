
############################################################################### -
# Classes ----------------------------------------------------------------------
############################################################################### -

#' An S4 class to store information
#'
#' @slot type (\code{character}) One of two available types "dense" or "sparse". See description.
#' @slot data (\code{data.frame}) A data.frame with covariates information, where colnames(data) matches model.list specification
#' @slot locs (\code{numeric matrix}) a matrix with locs matching data
#' @slot z (\code{numeric matrix}) A matrix of dimension n x p with response values
#' @slot model.list (\code{list}) A list specifying a model for each aspect of the spatial structure.
#' @slot info (\code{list}) a list with information about the coco object
#' @slot output (\code{list}) if building an already fitted \code{coco} object (not the standard approach), then requires an output from Optimparallel output, including as well boundaries, etc.
#' @author Federico Blasi
#' 
setClass("coco", slots = list(
  type = "character",
  data = "data.frame",
  locs = "matrix",
  z = "matrix",
  model.list = "list",
  info = "list",
  output = "list"
), package = "cocons", prototype = list(output = list()))

###############################################################################-
# Methods ----------------------------------------------------------------------
###############################################################################-

#' Plot Method for coco objects
#'
#' This method plots objects of class \code{coco}.
#' @rdname plot-methods
#' @docType methods
#' @aliases plot,coco-method
#' @param x (\code{S4}) A fitted object of class \code{coco}.
#' @param y Not used.
#' @param type (\code{character}  or \code{NULL}) The type of plot. NULL or "ellipse" for drawing ellipse of the convolution kernels.
#' @param index (\code{integer vector}) For plotting local correlation plots.
#' @param factr (\code{numeric}) Factor rate for size of ellipses.
#' @param ... Additional arguments passed to \link[fields]{quilt.plot}. 
#' @return Several plots are created.
#' @author Federico Blasi
#' 
setMethod("plot",
          signature(x = "coco", y = "missing"),
          definition =
            function(x, y, type = NULL, index = NULL, factr = 0.1, ...) {
              
              if (length(x@output) == 0) {
                stop("object has not yet been fitted.")
              }
              
              opar <- graphics::par(no.readonly = TRUE)
              on.exit(graphics::par(opar))
              
              if (any(x@type == "dense")) {

                spat_effects <- getSpatEffects(x)

                if(is.null(type)){
                  
                  graphics::par(mfrow = c(1, 2))
                  
                  do.call(fields::quilt.plot, list(
                    "x" = x@locs[, 1],
                    "y" = x@locs[, 2],
                    "main" = "spatial mean",
                    "z" = getSpatMean(x)
                  ), ...)
                  
                  if(T){
                    do.call(fields::quilt.plot, list(
                      "x" = x@locs[, 1],
                      "y" = x@locs[, 2],
                      "main" = "residuals",
                      "z" = x@z[,1] - getSpatMean(x)
                    ), ...)                    
                  }

                  graphics::par(mfrow = c(2, 3))
                  
                  fields::quilt.plot(x@locs, spat_effects$sd, main = "se", ...)
                  
                  #plot(x@locs[,1], x@locs[,2], pch=20, main = "se",
                  #     col = coco.pallete(64)[cut(spat_effects$sd, breaks = quantile(x = spat_effects$sd, probs = seq(0, 1, length.out = 64)), 
                  #                              labels = FALSE, include.lowest = T)], xlab = colnames(x@locs)[1],ylab=colnames(x@locs)[2])
                  
                  fields::quilt.plot(x@locs, spat_effects$scale_x, main = "approx. eff. scale", ...)
                  
                  #plot(x@locs[,1], x@locs[,2], pch=20, main = "approx eff.scale",
                  #     col = coco.pallete(64)[cut(spat_effects$scale_x, breaks = quantile(x = spat_effects$scale_x, probs = seq(0, 1, length.out = 64)), 
                  #                              labels = FALSE, include.lowest = T)], xlab = colnames(x@locs)[1],ylab=colnames(x@locs)[2])
                  
                  fields::quilt.plot(x@locs, spat_effects$aniso, main = "anisotropy", ...)
                  
                  #plot(x@locs[,1], x@locs[,2], pch=20, main = "ga",
                  #     col = coco.pallete(64)[cut(spat_effects$aniso, breaks = quantile(x = spat_effects$aniso, probs = seq(0, 1, length.out = 64)), 
                  #                              labels = FALSE, include.lowest = T)], xlab = colnames(x@locs)[1],ylab=colnames(x@locs)[2])
                  
                  fields::quilt.plot(x@locs, spat_effects$angle, main = "tilt", zlim = c(0, pi), ...)
                  fields::quilt.plot(x@locs, spat_effects$smooth, main = "smooth", ...)
                  fields::quilt.plot(x@locs, spat_effects$nugget, main = "nugget", ...)
                  
                }
                
                # ellipse grid
                
                if (!is.null(type)) {
                  if(type == "ellipse"){
                    
                    graphics::par(mfrow = c(1, 1))
                    
                    plot(x@locs, col = fields::tim.colors(128)[cut(x@z, 128)], pch = 20, cex = 1, asp = 1)
                    
                    number_x <- 10
                    range_x <- range(x@locs[, 1])
                    eps_x <- diff(range_x) / 3
                    vals_x <- seq(from = range_x[1], to = range_x[2], length.out = number_x)
                    
                    for (ii in 1:length(vals_x)) {
                      ref_y <- which(x@locs[, 1] < (vals_x[ii] + eps_x / 2) &
                                       x@locs[, 1] > (vals_x[ii] - eps_x / 2))
                      
                      range_y <- range(x@locs[ref_y, 2])
                      
                      for (jj in 1:number_x) {
                        center_locs <- c(vals_x[ii], range_y[1] + jj / number_x * (range_y[2] - range_y[1]))
                        
                        sss <- spam::nearest.dist(x = matrix(center_locs, ncol = 2), y = x@locs, delta = diff(range(x@locs))/5)
                        
                        if(length(sss@entries) == 0){next}
                        
                        to_compute_sd <- sss@colindices[which.min(sss@entries)]
                        
                        .cocons.DrawEllipsoid(
                          alpha_i = spat_effects$tilt[to_compute_sd], 
                          r = spat_effects$aniso[to_compute_sd], 
                          rho = spat_effects$scale_x[to_compute_sd],
                          loc = x@locs[to_compute_sd, ], factr = factr
                        )
                      }
                    }
                  }
                }
                
                if (!is.null(type)) {
                if (type == "correlations") {
                  
                  tmp_info <- cocons::getDesignMatrix(model.list = x@model.list, data = x@data)
                  
                  theta_list <- cocons::getModelLists(theta = x@output$par, 
                                                      par.pos = tmp_info$par.pos, type = "diff")
                  
                  X_std <- cocons::getScale(tmp_info$model.matrix,
                                            mean.vector = x@info$mean.vector,
                                            sd.vector = x@info$sd.vector
                  )
                  
                  tmp_cov <- stats::cov2cor(cocons::cov_rns(theta_list, 
                                                            locs = x@locs, 
                                                            x_covariates = X_std$std.covs,
                                                            smooth_limits = x@info$smooth.limits
                  ))
                  
                  graphics::par(mfrow = c(1, 2))
                  
                  for (ww in index) {
                    
                    fields::quilt.plot(x@locs, tmp_cov[ww, ], zlim = c(0, 1), main = paste0('global corr. at index ', ww))
                    graphics::points(x@locs[ww, 1], x@locs[ww, 2], pch = "X", bg = 'red', col = 'blue', cex = 2)
                    
                    local_var <- X_std$std.covs[ww, ] %*% theta_list$std.dev
                    local_mr <- X_std$std.covs[ww, ] %*% theta_list$scale
                    local_ga <- X_std$std.covs[ww, ] %*% theta_list$aniso
                    local_tt <- X_std$std.covs[ww, ] %*% theta_list$tilt
                    local_smtns <- diff(x@info$smooth.limits)/ (1 + exp(-X_std$std.covs[ww, ] %*% theta_list$smooth)) + x@info$smooth.limits[1]
                    local_nugget <- X_std$std.covs[ww, ] %*% theta_list$nugget
                    
                    model.list <- list(
                      "mean" = 0,
                      "std.dev" = stats::as.formula(" ~ 1"),
                      "scale" = stats::as.formula(" ~ 1"),
                      "aniso" = stats::as.formula(" ~ 1"),
                      "tilt" = stats::as.formula(" ~ 1"),
                      "smooth" = local_smtns,
                      "nugget" = -Inf
                    )
                    
                    local_object <- coco(
                      type = "dense",
                      locs = x@locs,
                      data = x@data,
                      z = x@z,
                      model.list = model.list
                    )
                    
                    DM_tmp <- getDesignMatrix(model.list = model.list, 
                                              data = local_object@data)
                    
                    here <- getModelLists(c(local_var, 
                                            local_mr, 
                                            local_ga, 
                                            local_tt, 
                                            local_smtns),
                                          par.pos = DM_tmp$par.pos, type = "diff"
                    )
                    
                    tmp_cov_two <- stats::cov2cor(cocons::cov_rns(
                      theta = here, locs = x@locs,
                      x_covariates = DM_tmp$model.matrix,
                      smooth_limits = x@info$smooth.limits
                    ))
                    
                    fields::quilt.plot(x@locs, tmp_cov_two[ww, ], main = paste0('local corr. at index ',ww))
                    graphics::points(x@locs[ww,1], x@locs[ww,2], pch = "X", bg='red', col = 'blue', cex = 2)
                  }
                }
                }
                
              }

              if (x@type == "sparse") {
                
                spat_effects <- getSpatEffects(x)
                
                tmp_list <- list("x" = x@locs[, 1],
                                 "y" = x@locs[, 2],
                                 "main" = "spatial mean",
                                 "z" = getSpatMean(x))
                
                graphics::par(mfrow = c(1, 2))
                do.call(fields::quilt.plot, args = tmp_list)
                
                tmp_list <- list("x" = x@locs[, 1],
                                 "y" = x@locs[, 2],
                                 "main" = "residuals",
                                 "z" = x@z[,1] - getSpatMean(x))
                
                do.call(fields::quilt.plot, args = tmp_list)
                
                graphics::par(mfrow = c(2, 2))
                
                fields::quilt.plot(x@locs, spat_effects$smooth, ..., main = "smooth")
                fields::quilt.plot(x@locs, spat_effects$scale_x, ..., main = "approx. eff. scale")
                fields::quilt.plot(x@locs, spat_effects$sd, ..., main = "se")
                fields::quilt.plot(x@locs, spat_effects$nugget, ..., main = "nugget")
                
                if (!is.null(type)) {
                  if(type == "ellipse"){
                    
                    graphics::par(mfrow = c(1, 1))
                    
                    plot(x@locs, col = fields::tim.colors(128)[cut(x@z, 128)], pch = 20, cex = 1, asp = 1)
                    
                    number_x <- 10
                    range_x <- range(x@locs[, 1])
                    eps_x <- diff(range_x) / 3
                    vals_x <- seq(from = range_x[1], to = range_x[2], length.out = number_x)
                    
                    for (ii in 1:length(vals_x)) {
                      ref_y <- which(x@locs[, 1] < (vals_x[ii] + eps_x / 2) &
                                       x@locs[, 1] > (vals_x[ii] - eps_x / 2))
                      
                      range_y <- range(x@locs[ref_y, 2])
                      
                      for (jj in 1:number_x) {
                        center_locs <- c(vals_x[ii], range_y[1] + jj / number_x * (range_y[2] - range_y[1]))
                        
                        sss <- spam::nearest.dist(x = matrix(center_locs, ncol = 2), y = x@locs, delta = diff(range(x@locs))/5) 
                        
                        to_compute_sd <- sss@colindices[which.min(sss@entries)]
                        
                        .cocons.DrawEllipsoid(
                          alpha_i = pi/2, 
                          r = 1, 
                          rho = spat_effects$scale_x[to_compute_sd],
                          loc = x@locs[to_compute_sd, ], factr = factr
                        )
                      }
                    }
                  }
                }
                
                if (!is.null(type)) {
                  if (type == "correlations") {
                  }
                }                
                
              }
            }
)

#' Summary Method for Coco Class
#'
#' method summary for objects of class 'coco'.
#' @name summary
#' @aliases summary,coco-method
#' @param object (\code{S4}) An object of class 'coco'.
#' @param inv.hess (\code{numeric matrix} or \code{NULL}) inverse of the approximated hessian matrix (getHessian)
#' @return summary the coco object
#' @docType methods
#' @rdname summary-methods
#' @author Federico Blasi
#' 
setMethod("summary", signature(object = "coco"), 
          definition = 
            function(object, inv.hess = NULL){
              
              if(length(object@output) == 0){stop("object has not been fited yet.")}
              
              tmp_matrix <- cocons::getDesignMatrix(model.list = object@model.list, 
                                                   data = object@data)
              
              adjusted_effects <- matrix(nrow = 7, ncol = dim(tmp_matrix$model.matrix)[2])
              colnames(adjusted_effects) <- colnames(tmp_matrix$model.matrix)
              
              adjusted_eff_values <- cocons::getModelLists(object@output$par, 
                                                          tmp_matrix$par.pos,
                                                          type = "diff") 

              if(!is.null(inv.hess)){
                
                Hess_mod <- getModHess(object, inv.hess = inv.hess)
                
                if(is.logical(tmp_matrix$par.pos$mean)){
                  number_mean <- sum(tmp_matrix$par.pos$mean)
                } else{number_mean <- 0}
                
                if(is.logical(tmp_matrix$par.pos$std.dev)){
                  number_std.dev <- sum(tmp_matrix$par.pos$std.dev)
                } else{number_std.dev <- 0}
                
                if(is.logical(tmp_matrix$par.pos$scale)){
                  number_scale <- sum(tmp_matrix$par.pos$scale)
                } else{number_scale <- 0}
                
                adjusted_se <- matrix(0,nrow = 7, ncol = dim(tmp_matrix$model.matrix)[2])
                colnames(adjusted_se) <- colnames(tmp_matrix$model.matrix)
                
                SE_diag <- sqrt(diag((Hess_mod)))
                cum_ses <- 0
                if(is.logical(tmp_matrix$par.pos$mean)){
                adjusted_se[1,tmp_matrix$par.pos$mean] <- SE_diag[1:number_mean]
                cum_ses <- number_mean
                }
                if(is.logical(tmp_matrix$par.pos$std.dev)){
                adjusted_se[2,tmp_matrix$par.pos$std.dev] <- SE_diag[(cum_ses+1):(cum_ses + number_std.dev)]
                cum_ses <- cum_ses + number_std.dev
                }
                if(is.logical(tmp_matrix$par.pos$scale)){
                adjusted_se[3,tmp_matrix$par.pos$scale] <- SE_diag[(cum_ses+1):(cum_ses + number_scale)]
                cum_ses <- cum_ses + number_scale
                }
                if(is.logical(tmp_matrix$par.pos$aniso)){
                adjusted_se[4,tmp_matrix$par.pos$aniso] <- SE_diag[(cum_ses+1):(cum_ses + sum(tmp_matrix$par.pos$aniso))]
                cum_ses <- cum_ses + sum(tmp_matrix$par.pos$aniso)
                }
                if(is.logical(tmp_matrix$par.pos$tilt)){
                adjusted_se[5,tmp_matrix$par.pos$tilt] <- SE_diag[(cum_ses+1):(cum_ses + sum(tmp_matrix$par.pos$tilt))]
                cum_ses <- cum_ses + sum(tmp_matrix$par.pos$tilt)
                }
                if(is.logical(tmp_matrix$par.pos$smooth)){
                  adjusted_se[6,tmp_matrix$par.pos$smooth] <- SE_diag[(cum_ses+1):(cum_ses + sum(tmp_matrix$par.pos$smooth))]
                  cum_ses <- cum_ses + sum(tmp_matrix$par.pos$smooth)
                }
                if(is.logical(tmp_matrix$par.pos$nugget)){
                  adjusted_se[7,tmp_matrix$par.pos$nugget] <- SE_diag[(cum_ses+1):(cum_ses + sum(tmp_matrix$par.pos$nugget))]
                }
                
              }
              
              if(object@type == "dense"){
                
                adjusted_effects[1, ] <- adjusted_eff_values$mean
                adjusted_effects[2, ] <- adjusted_eff_values$std.dev
                adjusted_effects[3, ] <- adjusted_eff_values$scale
                adjusted_effects[4, ] <- adjusted_eff_values$aniso
                adjusted_effects[5, ] <- adjusted_eff_values$tilt
                adjusted_effects[6, ] <- adjusted_eff_values$smooth
                adjusted_effects[7, ] <- adjusted_eff_values$nugget
                
                if(!is.null(inv.hess)){
                  
                  cat(sprintf("%-15s %15s %15s %15s\n", "Output:", " ", " ", " "))
                  cat(rep("-", 65), "\n")
                  
                  cat(sprintf("%-15s %15s %15s %15s %15s %15s %15s %15s\n", "(raw)", "Mean", "Std. Dev.", "Scale", "Aniso.", "Tilt", "Smooth", "Nugget"))
                  cat(rep("-", 65), "\n")
                  
                  for (ii in 1:dim(adjusted_effects)[2]) {
                    cat(sprintf("%-15s %15s %15s %15s %15s %15s %15s %15s\n", 
                                colnames(adjusted_effects)[ii], 
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$mean[ii]),
                                              yes = tmp_matrix$par.pos$mean[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                       paste0(round(adjusted_effects[1,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$mean[ii]), yes = round(adjusted_se[1, ii], 3),
                                                                                           no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$std.dev[ii]),
                                         yes = tmp_matrix$par.pos$std.dev[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[2,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$std.dev[ii]), yes = round(adjusted_se[2, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$scale[ii]),
                                         yes = tmp_matrix$par.pos$scale[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[3,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$scale[ii]), yes = round(adjusted_se[3, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$aniso[ii]),
                                         yes = tmp_matrix$par.pos$aniso[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[4,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$aniso[ii]), yes = round(adjusted_se[4, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$tilt[ii]),
                                         yes = tmp_matrix$par.pos$tilt[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[5,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$tilt[ii]), yes = round(adjusted_se[5, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$smooth[ii]),
                                         yes = tmp_matrix$par.pos$smooth[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[6,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$smooth[ii]), yes = round(adjusted_se[6, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$nugget[ii]),
                                         yes = tmp_matrix$par.pos$nugget[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[7,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$nugget[ii]), yes = round(adjusted_se[7, ii], 3),
                                                                                      no = "f"),")"), no = "-")
                    ))
                  }
                  
                } else {
                  
                  cat(sprintf("%-15s %15s %15s %15s\n", "Output:", " ", " ", " "))
                  cat(rep("-", 65), "\n")
                  
                  cat(sprintf("%-15s %10s %10s %10s %10s %10s %10s %10s\n", "(raw)", "Mean", "Std. Dev.", "Scale", "Aniso.", "Tilt", "Smooth", "Nugget"))
                  cat(rep("-", 65), "\n")
                  
                  for (ii in 1:dim(adjusted_effects)[2]) {
                    cat(sprintf("%-15s %10s %10s %10s %10s %10s %10s %10s\n", 
                                colnames(adjusted_effects)[ii], 
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$mean[ii]),
                                              yes = tmp_matrix$par.pos$mean[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[1,ii], 3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$std.dev[ii]),
                                              yes = tmp_matrix$par.pos$std.dev[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[2,ii], 3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$scale[ii]),
                                              yes = tmp_matrix$par.pos$scale[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[3,ii],3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$aniso[ii]),
                                              yes = tmp_matrix$par.pos$aniso[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[4,ii],3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$tilt[ii]),
                                              yes = tmp_matrix$par.pos$tilt[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[5,ii],3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$smooth[ii]),
                                              yes = tmp_matrix$par.pos$smooth[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[6, ii],3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$nugget[ii]),
                                              yes = tmp_matrix$par.pos$nugget[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[7, ii],3), "-")
                    )
                    )
                  }
                  
                }
                
                cat(rep("-", 65), "\n")
                
              }
              
              if(object@type == "sparse"){
                
                adjusted_effects[1, ] <- adjusted_eff_values$mean
                adjusted_effects[2, ] <- adjusted_eff_values$std.dev
                adjusted_effects[3, ] <- adjusted_eff_values$scale
                adjusted_effects[6, ] <- adjusted_eff_values$smooth
                adjusted_effects[7, ] <- adjusted_eff_values$nugget
                
                if(!is.null(inv.hess)){
                  
                  cat(sprintf("%-15s %15s %15s %15s\n", "Output:", " ", " ", " "))
                  cat(rep("-", 50), "\n")
                  
                  cat(sprintf("%-15s %15s %15s %15s %15s %15s\n", "(raw)", "Mean", "Std. Dev.", "Scale", "Smooth", "Nugget"))
                  cat(rep("-", 50), "\n")
                  
                  for (ii in 1:dim(adjusted_effects)[2]) {
                    cat(sprintf("%-15s %15s %15s %15s %15s %15s\n", 
                                colnames(adjusted_effects)[ii], 
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$mean[ii]),
                                         yes = tmp_matrix$par.pos$mean[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[1,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$mean[ii]), yes = round(adjusted_se[1, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$std.dev[ii]),
                                         yes = tmp_matrix$par.pos$std.dev[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[2,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$std.dev[ii]), yes = round(adjusted_se[2, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$scale[ii]),
                                         yes = tmp_matrix$par.pos$scale[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[3,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$scale[ii]), yes = round(adjusted_se[3, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$smooth[ii]),
                                         yes = tmp_matrix$par.pos$smooth[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[6,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$smooth[ii]), yes = round(adjusted_se[6, ii], 3),
                                                                                      no = "f"),")"), no = "-"),
                                ifelse(
                                  ifelse(is.logical(tmp_matrix$par.pos$nugget[ii]),
                                         yes = tmp_matrix$par.pos$nugget[ii], 
                                         no = ifelse(ii == 1, yes = TRUE, no = FALSE)), 
                                  paste0(round(adjusted_effects[7,ii], 3)," (",ifelse(is.logical(tmp_matrix$par.pos$smooth[ii]), yes = round(adjusted_se[7, ii], 3),
                                                                                      no = "f"),")"), no = "-")
                    )
                    )
                  }
                  
                } else {
                  
                  cat(sprintf("%-15s %15s %15s %15s\n", "Output:", " ", " ", " "))
                  cat(rep("-", 50), "\n")
                  
                  cat(sprintf("%-15s %10s %10s %10s %10s %10s\n", "(raw)", "Mean", "Std. Dev.", "Scale", "Smooth", "Nugget"))
                  cat(rep("-", 50), "\n")
                  
                  for (ii in 1:dim(adjusted_effects)[2]) {
                    cat(sprintf("%-15s %10s %10s %10s %10s %10s\n", 
                                colnames(adjusted_effects)[ii], 
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$mean[ii]),
                                              yes = tmp_matrix$par.pos$mean[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[1,ii], 3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$std.dev[ii]),
                                              yes = tmp_matrix$par.pos$std.dev[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[2,ii], 3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$scale[ii]),
                                              yes = tmp_matrix$par.pos$scale[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[3,ii], 3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$smooth[ii]),
                                              yes = tmp_matrix$par.pos$smooth[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[6, ii], 3), "-"),
                                ifelse(ifelse(is.logical(tmp_matrix$par.pos$nugget[ii]),
                                              yes = tmp_matrix$par.pos$nugget[ii], 
                                              no = ifelse(ii == 1, yes = TRUE, no = FALSE)), round(adjusted_effects[7, ii], 3), "-")
                    )
                    )
                  }
                  
                }
                
                cat(rep("-", 50), "\n")
                
              
              }
            }
)



#' Show Method for Coco Class
#'
#' This method show objects of class 'coco'.
#' @name show
#' @aliases show,coco-method
#' @param object (\code{S4}) An object of class 'coco'.
#' @return A plot is created.
#' @docType methods
#' @rdname show-methods
#' @author Federico Blasi
#' 
setMethod("show",
          signature(object = "coco"),
          function(object){
            
            if(object@type == "dense"){
              
              cat(sprintf("%-30s %30s\n", "coco object", object@type))
              cat(sprintf("%-30s %30s\n", "fitted", ifelse(identical(object@output, list()), yes = "no", 
                                                           no = "yes")))
              
              cat(rep("-", 40), "\n")
              cat(sprintf("%-30s %30s\n", "dataset dim", paste0(dim(object@data), collapse = ", ")))
              
              # number of parameters
              
              cat(sprintf("%-30s %30s\n", "model parameters", .cocons.get.npars(object)))
              
              # covariates names
              
              coco_info_light <- getDesignMatrix(model.list = object@model.list, data = object@data[1, , drop = FALSE])
              
              cat(sprintf("%-30s %10s", "covariates ", paste0(colnames(coco_info_light$model.matrix),collapse = ", ")))
              
            }
            
            if(object@type == "sparse"){
              
              cat(sprintf("%-30s %30s\n", "coco object", object@type))
              cat(sprintf("%-30s %30s\n", "fitted", ifelse(identical(object@output, list()), yes = "no", 
                                                           no = "yes")))
              
              cat(rep("-", 40), "\n")
              cat(sprintf("%-30s %30s\n", "dataset dim", paste0(dim(object@data), collapse = ", ")))
              
              # number of parameters
              
              cat(sprintf("%-30s %30s\n", "model parameters", .cocons.get.npars(object)))
              
              # covariates names
              
              coco_info_light <- getDesignMatrix(model.list = object@model.list, data = object@data[1, , drop = FALSE])
              
              cat(sprintf("%-30s %10s", "covariates ", paste0(colnames(coco_info_light$model.matrix),collapse = ", ")))
              
            }
            
          }
)
