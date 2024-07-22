
############################################################################### -
# Classes ----------------------------------------------------------------------
############################################################################### -

#' An S4 class to store information
#'
#' @slot type One of two available types "dense" or "sparse". See description.
#' @slot data A data.frame with covariates information, where colnames(data) matches model.list specification
#' @slot locs a matrix with locs matching data
#' @slot z A vector with response values
#' @slot model.list A list specyfing a model for each aspect of the spatial structure.
#' @slot info a list with information about the coco object
#' @slot output an output from optimparallel output, including as well boundaries 
#' information as another element of the list
#' @author Federico Blasi
#' 
setClass("coco", slots = list(
  type = "character",
  data = "data.frame",
  locs = "matrix",
  z = "matrix",
  model.list = "list",
  info = 'list',
  output = "list"
), package = "coco", prototype = list(output = list()))

###############################################################################-
# Methods ----------------------------------------------------------------------
###############################################################################-

#' Plot Method for Coco Class
#'
#' This method plots objects of class \code{coco}.
#'
#' @param x An object of class \code{coco}.
#' @param y Not used.
#' @param ... Additional arguments passed to the plot function.
#' @param type The type of plot.
#' @param index For plotting local correlation plots.
#' @param factr Factor rate for size of ellipses.
#' @param plot.control Additional plot control parameters.
#' @return A plot is created.
#' @exportMethod plot
#' @docType methods
#' @rdname plot-methods
#' @aliases plot,coco-method
setMethod("plot",
          signature(x = "coco", y = "missing"),
          definition =
            function(x, y, ..., type = NULL, index = NULL, factr = 0.1, plot.control = NULL) {
              
              if (length(x@output) == 0) {
                stop("object has not yet been fitted.")
              }
              
              opar <- graphics::par(no.readonly = TRUE)
              on.exit(graphics::par(opar))
              
              if (any(x@type == "dense")) {
                tmp_info <- cocons::getDesignMatrix(model.list = x@model.list, 
                                                   data = x@data)
                
                theta_list <- cocons::getModelLists(
                  theta = x@output$par, par.pos = tmp_info$par.pos,
                  type = "diff"
                )
                
                X_std <- cocons::getScale(tmp_info$model.matrix,
                                         mean.vector = x@info$mean.vector,
                                         sd.vector = x@info$sd.vector
                )
                
                tp_se <- exp(0.5 * X_std$std.covs %*% theta_list$std.dev)
                tp_ga <- exp(X_std$std.covs %*% theta_list$aniso)
                tp_tl <- pi / (1 + exp(-X_std$std.covs %*% theta_list$tilt))
                tp_smooth <- (x@info$smooth_limits[2] - x@info$smooth_limits[1]) / (1 + exp(-X_std$std.covs %*% theta_list$smooth)) + x@info$smooth_limits[1]
                tp_ng <- exp(X_std$std.covs %*% theta_list$nugget)
                tp_mr <- sin(tp_tl) * exp(X_std$std.covs %*% theta_list$scale) # *  exp(X_std$std.covs %*% theta_list$aniso)
                
                if(is.null(type)){
                  
                  graphics::par(mfrow = c(1, 2))
                  
                  do.call(fields::quilt.plot, list(
                    "x" = x@locs[, 1],
                    "y" = x@locs[, 2],
                    "main" = "trend",
                    "z" = X_std$std.covs %*% theta_list$mean
                  ), ...)
                  
                  if(T){
                    do.call(fields::quilt.plot, list(
                      "x" = x@locs[, 1],
                      "y" = x@locs[, 2],
                      "main" = "residuals",
                      "z" = x@z[,1] - X_std$std.covs %*% theta_list$mean
                    ), ...)                    
                  }

                  graphics::par(mfrow = c(2, 3))
                  
                  fields::quilt.plot(x@locs, tp_se, main = "se", ...)
                  
                  #plot(x@locs[,1], x@locs[,2], pch=20, main = 'se',
                  #     col = tim.colors(64)[cut(tp_se, breaks = quantile(x = tp_se, probs = seq(0, 1, length.out = 64)), 
                  #                              labels = FALSE, include.lowest = T)], xlab = colnames(x@locs)[1],ylab=colnames(x@locs)[2])
                  
                  fields::quilt.plot(x@locs, tp_mr, main = "eff.range", ...)
                  
                  #plot(x@locs[,1], x@locs[,2], pch=20, main = 'eff.range',
                  #     col = tim.colors(64)[cut(tp_mr, breaks = quantile(x = tp_mr, probs = seq(0, 1, length.out = 64)), 
                  #                              labels = FALSE, include.lowest = T)], xlab = colnames(x@locs)[1],ylab=colnames(x@locs)[2])
                  
                  fields::quilt.plot(x@locs, tp_ga, main = "ga", ...)
                  
                  #plot(x@locs[,1], x@locs[,2], pch=20, main = 'ga',
                  #     col = tim.colors(64)[cut(tp_ga, breaks = quantile(x = tp_ga, probs = seq(0, 1, length.out = 64)), 
                  #                              labels = FALSE, include.lowest = T)], xlab = colnames(x@locs)[1],ylab=colnames(x@locs)[2])
                  
                  fields::quilt.plot(x@locs, cos(tp_tl), main = "tl", zlim = c(-pi / 4 - 0.1, pi / 4 + 0.1), ...)
                  fields::quilt.plot(x@locs, tp_smooth, main = "smooth", ...)
                  fields::quilt.plot(x@locs, tp_ng, main = "ng", ...)
                  
                }
                
                # ellipse grid
                
                if (!is.null(type)) {
                  if(type == 'ellipse'){
                    
                    graphics::par(mfrow = c(1, 1))
                    
                    plot(x@locs, col = fields::tim.colors(128)[cut(x@z, 128)], pch = 18, cex = 2, asp = 1)
                    
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
                        
                        sss <- spam::nearest.dist(x = matrix(center_locs, ncol = 2), y = x@locs, delta = 2) # fix delta to automatic
                        
                        to_compute_sd <- sss@colindices[which.min(sss@entries)]
                        
                        .coco.DrawEllipsoid(
                          alpha_i = tp_tl[to_compute_sd], r = tp_ga[to_compute_sd], rho = tp_mr[to_compute_sd],
                          loc = x@locs[to_compute_sd, ], factr = factr
                        )
                      }
                    }
                  }
                }
                
              }
              
              if (!is.null(type)) {
                
                # !!!!!!! is it ok && ??
                if (x@type == "dense" && type == "correlations") {

                  tmp_cov <- stats::cov2cor(cocons::cov_rns(theta_list, x@locs, X_std$std.covs,
                                                           smooth_limits = x@info$smooth_limits
                  ))
                  
                  theta_list <- cocons::getModelLists(
                    theta = x@output$par, par.pos = tmp_info$par.pos,
                    type = "classic"
                  )
                  
                  # ifelse(length(index) <= 5, yes = graphics::par(mfrow = c(length(index), 2)), no = graphics::par(mfrow = c(5, 2)))

                  graphics::par(mfrow = c(1, 2))
                  
                  for (ww in index) {
                    
                    fields::quilt.plot(x@locs, tmp_cov[ww,], zlim = c(0, 1))
                    local_var <- X_std$std.covs[ww, ] %*% theta_list$std.dev
                    local_mr <- X_std$std.covs[ww, ] %*% theta_list$scale
                    local_ga <- X_std$std.covs[ww, ] %*% theta_list$aniso
                    local_tt <- X_std$std.covs[ww, ] %*% theta_list$tilt
                    local_smtns <- X_std$std.covs[ww, ] %*% theta_list$smooth
                    local_nugget <- X_std$std.covs[ww, ] %*% theta_list$nugget
                    
                    model.list <- list(
                      "mean" = 0,
                      "std.dev" = stats::as.formula(" ~ 1"),
                      "scale" = stats::as.formula(" ~ 1"),
                      "aniso" = stats::as.formula(" ~ 1"),
                      "tilt" = stats::as.formula(" ~ 1"),
                      "smooth" = stats::as.formula(" ~ 1"),
                      "nugget" = -Inf
                    )
                    
                    local_object <- coco(
                      type = "dense",
                      locs = x@locs,
                      data = x@data,
                      z = x@z,
                      model.list = model.list,
                      info = list("smooth_limits" = x@info$smooth_limits)
                    )
                    
                    teteee <- getDesignMatrix(model.list = model.list, 
                                              data = local_object@data)
                    
                    here <- getModelLists(c(local_var, 
                                            local_mr, 
                                            local_ga, 
                                            local_tt, 
                                            local_smtns),
                                          par.pos = teteee$par.pos, type = "diff"
                    )
                    
                    tmp_cov_two <- stats::cov2cor(cocons::cov_rns(
                      theta = here, locs = x@locs,
                      x_covariates = teteee$model.matrix,
                      smooth_limits = x@info$smooth_limits
                    ))
                    
                    fields::quilt.plot(x@locs, tmp_cov_two[ww, ])
                  }
                }
              }
              
              if (x@type == "sparse") {
                
                tmp_info <- cocons::getDesignMatrix(model.list = x@model.list, data = x@data)
                
                theta_list <- cocons::getModelLists(theta = x@output$par, 
                                                   par.pos = tmp_info$par.pos, type = "diff")
                
                X_std <- cocons::getScale(tmp_info$model.matrix,
                                         mean.vector = x@info$mean.vector,
                                         sd.vector = x@info$sd.vector
                )
                
                tp_smooth <- (x@info$smooth_limits[2] - x@info$smooth_limits[1]) / (1 + exp(-X_std$std.covs %*% theta_list$smooth)) + x@info$smooth_limits[1]
                
                # tp_eff_range <- sqrt(8 * tp_smooth) * exp(X_std$std.covs %*% theta_list$scale)
                tp_eff_range <- exp(X_std$std.covs %*% theta_list$scale)
                
                tp_se <- exp(0.5 * X_std$std.covs %*% theta_list$std.dev)
                
                tp_nugget_se <- sqrt(exp(X_std$std.covs %*% theta_list$nugget))
                
                tmp_list <- list('x' = x@locs[, 1],
                                 'y' = x@locs[, 2],
                                 'main' = 'trend',
                                 'z' = X_std$std.covs %*% theta_list$mean)
                
                graphics::par(mfrow = c(1, 2))
                do.call(fields::quilt.plot, args = tmp_list)
                
                tmp_list <- list('x' = x@locs[, 1],
                                 'y' = x@locs[, 2],
                                 'main' = 'residuals',
                                 'z' = x@z - X_std$std.covs %*% theta_list$mean)
                
                do.call(fields::quilt.plot, args = tmp_list)
                
                graphics::par(mfrow = c(2, 2))
                
                fields::quilt.plot(x@locs, tp_smooth, ..., main = "smooth")
                fields::quilt.plot(x@locs, tp_eff_range, ..., main = "approx eff.range")
                fields::quilt.plot(x@locs, tp_se, ..., main = "se")
                fields::quilt.plot(x@locs, tp_nugget_se, ..., main = "nugget se")
              }
              
              if (!is.null(type)) {
                if (x@type == "sparse" & type == "correlations") {
                  tmp_info <- cocons::getDesignMatrix(
                    model.list = x@model.list,
                    data = x@data
                  )
                  
                  theta_list <- cocons::getModelLists(
                    theta = x@output$par,
                    par.pos = tmp_info$par.pos, type = "diff"
                  )
                  
                  X_std <- cocons::getScale(tmp_info$model.matrix)
                  
                  tmp_cov <- stats::cov2cor(cocons::cov_rns_smooth_taper_vector(
                    theta_list, x@locs,
                    X_std$std.covs
                  )) ## !!!!!!!!!!!!!!!!!!!!!!!!!! FIXXXXXXXX
                  
                  for (ww in 1:length(index)) {
                    fields::quilt.plot(x@locs, tmp_cov[index[ww], ])
                  }
                }
              }
            }
)

#' Print Method for Coco Class
#'
#' This method prints objects of class 'coco'.
#' @name print
#' @aliases print,coco-method
#' @param x An object of class 'coco'.
#' @param inv.hess inverse of the approximated hessian matrix (getHessian)
#' @param ... Additional arguments to be passed to plot.
#' @return print the coco object
#' @docType methods
#' @exportMethod print
#' @docType methods
#' @rdname print-methods
#' @author Federico Blasi
#' 
setMethod("print", signature(x = "coco"), 
          definition = 
            function(x, inv.hess = NULL, ...){
              
              if(length(x@output) == 0){stop('object has not been fited yet.')}
              
              tmp_matrix <- cocons::getDesignMatrix(model.list = x@model.list, 
                                                   data = x@data)
              
              adjusted_effects <- matrix(nrow = 7, ncol = dim(tmp_matrix$model.matrix)[2])
              colnames(adjusted_effects) <- colnames(tmp_matrix$model.matrix)
              
              adjusted_eff_values <- cocons::getModelLists(x@output$par, 
                                                          tmp_matrix$par.pos,
                                                          type = 'diff') 

              if(!is.null(inv.hess)){
                
                Hess_mod <- getModHess(x, inv.hess = inv.hess)
                
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
              
              if(x@type == 'dense'){
                
                adjusted_effects[1, ] <- adjusted_eff_values$mean
                adjusted_effects[2, ] <- adjusted_eff_values$std.dev
                adjusted_effects[3, ] <- adjusted_eff_values$scale
                adjusted_effects[4, ] <- adjusted_eff_values$aniso
                adjusted_effects[5, ] <- adjusted_eff_values$tilt
                adjusted_effects[6, ] <- adjusted_eff_values$smooth
                adjusted_effects[7, ] <- adjusted_eff_values$nugget
                

                #cat(sprintf("%-15s %-15s %-15s %-15s\n", "Variable:", " "," "," "))
                
                if(!is.null(inv.hess)){
                  
                  # Create a table-like output using cat() and sprintf()
                  cat(sprintf("%-15s %15s %15s %15s\n", "Output:", " ", " ", " "))
                  cat(rep("-", 65), "\n")
                  
                  cat(sprintf("%-15s %15s %15s %15s %15s %15s %15s %15s\n", "(raw)", "Mean", "Std. Dev.", "Scale", "Geom. Aniso.", "Tilt", "Smooth", "Nugget"))
                  cat(rep("-", 65), "\n")
                  
                  for (ii in 1:dim(adjusted_effects)[2]) {
                    cat(sprintf("%-15s %15s %15s %15s %15s %15s %15s %15s\n", 
                                colnames(adjusted_effects)[ii], 
                                ifelse(tmp_matrix$par.pos$mean[ii], 
                                       paste0(round(adjusted_effects[1,ii], 3)," (",round(adjusted_se[1,ii],3),")"), "-"),
                                ifelse(tmp_matrix$par.pos$std.dev[ii], yes = paste0(round(adjusted_effects[2,ii], 3)," (",round(adjusted_se[2,ii], 3),")"), 
                                       no = "-"),
                                ifelse(tmp_matrix$par.pos$scale[ii], yes = paste0(round(adjusted_effects[3,ii],3)," (",round(adjusted_se[3,ii], 3),")"), 
                                       no = "-"),
                                ifelse(tmp_matrix$par.pos$aniso[ii], yes = paste0(round(adjusted_effects[4,ii],3)," (",round(adjusted_se[4,ii], 3),")"), 
                                       no = "-"),
                                ifelse(tmp_matrix$par.pos$tilt[ii], yes = paste0(round(adjusted_effects[5,ii],3)," (",round(adjusted_se[5,ii], 3),")"), 
                                       no = "-"),
                                ifelse(tmp_matrix$par.pos$smooth[ii], yes = paste0(round(adjusted_effects[6, ii],3)," (",round(adjusted_se[6,ii], 3),")"), 
                                       no = "-"),
                                ifelse(!is.na(tmp_matrix$par.pos$nugget[ii]), ifelse(tmp_matrix$par.pos$nugget[ii], paste0(round(adjusted_effects[7,ii], 3)," (",
                                                                                                                           round(adjusted_se[7,ii], 3),")"), "-"),
                                                                                     no = '-')
                    )
                    )
                  }
                  
                } else {
                  
                  # Create a table-like output using cat() and sprintf()
                  cat(sprintf("%-15s %15s %15s %15s\n", "Output:", " ", " ", " "))
                  cat(rep("-", 60), "\n")
                  
                  cat(sprintf("%-15s %10s %10s %10s %10s %10s %10s %10s\n", "(raw)", "Mean", "Std. Dev.", "Scale", "Geom. Aniso.", "Tilt", "Smooth", "Nugget"))
                  cat(rep("-", 60), "\n")
                  
                  for (ii in 1:dim(adjusted_effects)[2]) {
                    cat(sprintf("%-15s %10s %10s %10s %10s %10s %10s %10s\n", 
                                colnames(adjusted_effects)[ii], 
                                ifelse(tmp_matrix$par.pos$mean[ii], round(adjusted_effects[1,ii], 3), "-"),
                                ifelse(tmp_matrix$par.pos$std.dev[ii], round(adjusted_effects[2,ii], 3), "-"),
                                ifelse(tmp_matrix$par.pos$scale[ii], round(adjusted_effects[3,ii],3), "-"),
                                ifelse(tmp_matrix$par.pos$aniso[ii], round(adjusted_effects[4,ii],3), "-"),
                                ifelse(tmp_matrix$par.pos$tilt[ii], round(adjusted_effects[5,ii],3), "-"),
                                ifelse(tmp_matrix$par.pos$smooth[ii], round(adjusted_effects[6, ii],3), "-"),
                                ifelse(!is.na(tmp_matrix$par.pos$nugget[ii]), ifelse(tmp_matrix$par.pos$nugget[ii], round(adjusted_effects[7,ii],3),'-'),'-')
                    )
                    )
                  }
                  
                }
                
                cat(rep("-", 65), "\n")
                
              }
              
              if(x@type == "sparse"){
                
                # warning('taper estimations are biased.')
                
                adjusted_effects[1, ] <- adjusted_eff_values$mean
                adjusted_effects[2, ] <- adjusted_eff_values$std.dev
                adjusted_effects[3, ] <- adjusted_eff_values$scale
                adjusted_effects[6, ] <- adjusted_eff_values$smooth
                adjusted_effects[7, ] <- adjusted_eff_values$nugget
                
                
                #cat(sprintf("%-15s %-15s %-15s %-15s\n", "Variable:", " "," "," "))
                
                if(!is.null(inv.hess)){
                  
                  # Create a table-like output using cat() and sprintf()
                  cat(sprintf("%-15s %15s %15s %15s\n", "Output:", " ", " ", " "))
                  cat(rep("-", 65), "\n")
                  
                  cat(sprintf("%-15s %15s %15s %15s %15s %15s\n", "(raw)", "Mean", "Std. Dev.", "Scale", "Smooth", "Nugget"))
                  cat(rep("-", 65), "\n")
                  
                  
                  for (ii in 1:dim(adjusted_effects)[2]) {
                    cat(sprintf("%-15s %15s %15s %15s %15s %15s\n", 
                                colnames(adjusted_effects)[ii], 
                                ifelse(tmp_matrix$par.pos$mean[ii], 
                                       paste0(round(adjusted_effects[1,ii], 3)," (",round(adjusted_se[1,ii],3),")"), "-"),
                                ifelse(tmp_matrix$par.pos$std.dev[ii], yes = paste0(round(adjusted_effects[2,ii], 3)," (",round(adjusted_se[2,ii], 3),")"), 
                                       no = "-"),
                                ifelse(tmp_matrix$par.pos$scale[ii], yes = paste0(round(adjusted_effects[3,ii],3)," (",round(adjusted_se[3,ii], 3),")"), 
                                       no = "-"),
                                ifelse(tmp_matrix$par.pos$smooth[ii], yes = paste0(round(adjusted_effects[6, ii],3)," (",round(adjusted_se[6,ii], 3),")"), 
                                       no = "-"),
                                ifelse(!is.na(tmp_matrix$par.pos$nugget[ii]), ifelse(tmp_matrix$par.pos$nugget[ii], paste0(round(adjusted_effects[7,ii], 3)," (",
                                                                                                                           round(adjusted_se[7,ii], 3),")"), "-"),
                                       no = '-')
                    )
                    )
                  }
                  
                } else {
                  
                  # Create a table-like output using cat() and sprintf()
                  cat(sprintf("%-15s %15s %15s %15s\n", "Output:", " ", " ", " "))
                  cat(rep("-", 60), "\n")
                  
                  cat(sprintf("%-15s %10s %10s %10s %10s %10s\n", "(raw)", "Mean", "Std. Dev.", "Scale", "Smooth", "Nugget"))
                  cat(rep("-", 60), "\n")
                  
                  for (ii in 1:dim(adjusted_effects)[2]) {
                    cat(sprintf("%-15s %10s %10s %10s %10s %10s\n", 
                                colnames(adjusted_effects)[ii], 
                                ifelse(tmp_matrix$par.pos$mean[ii], round(adjusted_effects[1,ii], 3), "-"),
                                ifelse(tmp_matrix$par.pos$std.dev[ii], round(adjusted_effects[2,ii], 3), "-"),
                                ifelse(tmp_matrix$par.pos$scale[ii], round(adjusted_effects[3,ii],3), "-"),
                                ifelse(tmp_matrix$par.pos$smooth[ii], round(adjusted_effects[6, ii],3), "-"),
                                ifelse(!is.na(tmp_matrix$par.pos$nugget[ii]), ifelse(tmp_matrix$par.pos$nugget[ii], round(adjusted_effects[7,ii],3),'-'),'-')
                    )
                    )
                  }
                  
                }
                
                cat(rep("-", 65), "\n")
                
              
              }
            }
)

#' Show Method for Coco Class
#'
#' This method show objects of class 'coco'.
#' @name show
#' @aliases show,coco-method
#' @param object An object of class 'coco'.
#' @return A plot is created.
#' @docType methods
#' @rdname show-methods
#' @exportMethod show
#' @author Federico Blasi
#' 
setMethod("show",
          signature(object = "coco"),
          function(object){
            
            if(object@type == 'dense'){
              cat(sprintf("%-30s %30s\n", "coco object", object@type))
              cat(sprintf("%-30s %30s\n", "fitted", ifelse(identical(object@output, list()), yes = 'no', 
                                                           no = 'yes')))
              
              cat(rep("-", 40), "\n")
              cat(sprintf("%-30s %30s\n", "dataset dim", paste0(dim(object@data), collapse = ', ')))
              
              # number of parameters
              
              coco_info_light <- getDesignMatrix(model.list = object@model.list,data = object@data[1:2, ])
              
              tmp_index <- lapply(coco_info_light$par.pos, FUN = is.logical)
              
              n_par <- sum(unlist(lapply(coco_info_light$par.pos, sum))[which(tmp_index == TRUE)])
              
              cat(sprintf("%-30s %30s\n", "model parameters", n_par))
              
              # covariates names
              
              cat(sprintf("%-30s %10s", "covariates ", paste0(colnames(coco_info_light$model.matrix),collapse = ', ')))
              
            }
            
            if(object@type == 'sparse'){
              cat(sprintf("%-30s %30s\n", "coco object", object@type))
              cat(sprintf("%-30s %30s\n", "fitted", ifelse(identical(object@output, list()), yes = 'no', 
                                                           no = 'yes')))
              # cat(sprintf("%-30s %30s\n", "taper function", paste(object@info$taper)))
              
              cat(rep("-", 40), "\n")
              cat(sprintf("%-30s %30s\n", "dataset dim", paste0(dim(object@data), collapse = ', ')))
              
              # number of parameters
              
              coco_info_light <- getDesignMatrix(model.list = object@model.list,data = object@data[1:2, ])
              
              tmp_index <- lapply(coco_info_light$par.pos, FUN = is.logical)
              
              n_par <- sum(unlist(lapply(coco_info_light$par.pos, sum))[which(tmp_index == TRUE)])
              
              cat(sprintf("%-30s %30s\n", "model parameters", n_par))
              
              # covariates names
              
              cat(sprintf("%-30s %10s", "covariates ", paste0(colnames(coco_info_light$model.matrix),collapse = ', ')))
              
            }
            
          }
)
