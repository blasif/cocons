
# Internal functions

.cocons.update.optim.control <- function(optim.control){
  
  to_update <- getOption("cocons.Optim.Control")
  
  if(is.null(optim.control$parallel$forward)){
    optim.control$parallel$forward <- to_update$parallel$forward
  }
  
  if(is.null(optim.control$parallel$loginfo)){
    optim.control$parallel$loginfo <- to_update$parallel$loginfo
  }
  
  if(is.null(optim.control$control$factr)){
    optim.control$control$factr <- to_update$control$factr
  }
  
  if(is.null(optim.control$control$ndeps)){
    optim.control$control$ndeps <- to_update$control$ndeps
  }
  
  if(is.null(optim.control$control$trace)){
    optim.control$control$trace <- to_update$control$trace
  }
  
  if(is.null(optim.control$control$maxit)){
    optim.control$control$maxit <- to_update$control$maxit
  }
  
  if(is.null(optim.control$control$lmm)){
    optim.control$control$lmm <- to_update$control$lmm
  }
  
  optim.control$hessian <- FALSE

  return(optim.control)
  
}

.cocons.DrawEllipsoid <- function(alpha_i, r, rho, loc, factr) {
  
  createEllipse <- function(center, a, b, angle = 0, steps = 100, factr = 0.1) {
    theta <- seq(0, 2 * pi, length.out = steps)
    ellipse_points <- cbind(a * cos(theta), b * sin(theta))
    ellipse_points <- ellipse_points %*% matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
    ellipse_points[, 1] <- factr * ellipse_points[, 1] + center[1]
    ellipse_points[, 2] <- factr * ellipse_points[, 2] + center[2]
    
    return(ellipse_points)
  }
  
  if (abs(alpha_i - pi / 2) <= 1e-3) {
    
    if(r == 1){
      
      if(rho * factr > 0.05){
        
        graphics::arrows(
          x0 = loc[1], x1 = loc[1] + rho * factr,
          y0 = loc[2], y1 = loc[2], lwd = 2, lty = 1,
          cex = 0.5, angle = 5, length = 0.1
        )
        
        graphics::arrows(
          x0 = loc[1], x1 = loc[1],
          y0 = loc[2], y1 = loc[2] + rho * factr, lwd = 2, lty = 1,
          cex = 0.5, angle = 5, length = 0.1
        )
      }
      
      ellipse_points <- createEllipse(loc,
                                      a = rho,
                                      b = rho, 0, steps = 2000,
                                      factr = factr
      )
      
      graphics::lines(ellipse_points, lty = 3, col = "black")      
      
      return(0)
    }
    
    if(rho * factr > 0.05){
     
      graphics::arrows(
        x0 = loc[1], x1 = loc[1] + rho * factr,
        y0 = loc[2], y1 = loc[2] , lwd = 2, lty = 1,
        cex = 0.5, angle = 5, length = 0.1
      )
      
      graphics::arrows(
        x0 = loc[1], x1 = loc[1] ,
        y0 = loc[2], y1 = loc[2] + rho * r * factr, lwd = 2, lty = 1,
        cex = 0.5, angle = 5, length = 0.1
      )      
       
    }
    
    ellipse_points <- createEllipse(loc,
                                    a = max(rho,rho * r),
                                    b = min(rho,rho * r), 0, steps = 2000,
                                    factr = factr
    )
    
    graphics::lines(ellipse_points, lty = 3, col = "black")
    return(0)
    
  }
  
  A_stuff <- sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2)
  
  e_1 <- rho / 2 * ((r + 1) + sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
  e_2 <- rho / 2 * ((r + 1) - sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
  
  a_a_1 <- 2 * r^(0.5) * cospi(alpha_i / pi)
  a_b_1 <- r - 1 - A_stuff
  
  e_b_1 <- 2 * r^(0.5) * cospi(alpha_i / pi)
  e_b_2 <- r - 1 + A_stuff
  
  magnitude_vector_one <- sqrt((a_a_1)^2 + (a_b_1)^2)
  
  y_length <- factr * sqrt(e_2) / magnitude_vector_one
  
  magnitude_vector_two <- sqrt((e_b_1)^2 + (e_b_2)^2)
  
  y_length_2 <- factr * sqrt(e_1) / magnitude_vector_two
  
  if( (y_length * a_a_1 > 0.05) & (y_length_2 * e_b_1 > 0.05)){

      graphics::arrows(
    x0 = loc[1], x1 = loc[1] + y_length * a_a_1,
    y0 = loc[2], y1 = loc[2] + y_length * a_b_1, lwd = 2, lty = 1,
    cex = 0.5, angle = 15, length = 0.1, col = "black"
    )
    
    graphics::arrows(
      x0 = loc[1], x1 = loc[1] + y_length_2 * e_b_1,
      y0 = loc[2], y1 = loc[2] + y_length_2 * e_b_2, lwd = 2, lty = 1,
      cex = 0.5, angle = 15, length = 0.1, col = "black"
    )
    
  }
  
  ellipse_points <- createEllipse(
    center = loc, a = sqrt((max(e_1, e_2))),
    b = sqrt((min(e_1, e_2))), atan2(x = e_b_1, y = e_b_2),
    steps = 2000, factr = factr
  )
  
  graphics::lines(ellipse_points, lty = 3, col = "black")
}

.cocons.check.type_pred <- function(type){
  if(!(type %in% c("mean", "pred"))){
    stop(" type must be `mean` or `pred`")
  }
}

.cocons.get.npars <- function(coco.object){
  
  coco_info_light <- getDesignMatrix(model.list = coco.object@model.list, data = coco.object@data[1, , drop = FALSE])
  
  tmp_index <- lapply(coco_info_light$par.pos, FUN = is.logical)
  
  n_par <- sum(unlist(lapply(coco_info_light$par.pos, sum))[which(tmp_index == TRUE)])
  
  return(n_par)

}

.cocons.check.pars <- function(coco.object, pars){
  
  if(is.null(pars) & length(coco.object@output) > 1){return(0)}
  
  stopifnot("length of pars does not match number of parameters to estimate in the coco object." =
              .cocons.get.npars(coco.object) == length(pars))
}

.cocons.set.ncores <- function(coco.object, optim.control){
  
  n_pars <- .cocons.get.npars(coco.object)
  n_threads <- parallel::detectCores()
  
  if(n_threads < 3){
    return( n_threads)
  }

  if(is.null(optim.control)){
    n_total <- 2 * n_pars + 1
  } else{
    if(!is.null(optim.control$parallel$forward)){
      if(optim.control$parallel$forward){
        n_total <- n_pars + 1
      }else{
      n_total <- 2 * n_pars + 1
    }
    } else {n_total <- 2 * n_pars + 1}
  }
  
  n_total_t <- n_total
  p <- 1
  while(T){
    if(n_threads > n_total_t){return(n_total_t)}else{
      p <- p + 1
      n_total_t <- ceiling(n_total/p)
    }
  }
}

.cocons.check.type <- function(type) {
  if (!(type %in% c("sparse", "dense"))) {
    stop("check type")
  }
  return(0)
}

.cocons.check.coco <- function(coco){
  stopifnot("not a coco object" = methods::is(coco, 'coco'))
}

.cocons.check.data <- function(data) {

  stopifnot("data must be provided as data.frame" = is.data.frame(data),
            "data does not have colnames, which are needed for building `par.pos` " = !is.null(colnames(data)))

  return(0)
}

.cocons.check.locs <- function(locs) {
  
  stopifnot("locs should be provided as matrix" = is.matrix(locs))

  return(0)
}

.cocons.check.z <- function(z, data) {
  if (is.null(z)) {
    warning("z not provided. Expecting to simulate with this coco object.")
  } else {
    stopifnot("z is not a matrix." =  is.matrix(z),
              "dim(z)[1] != dim(data)[1]" = dim(z)[1] == dim(data)[1])
  }
  
  return(0)
}

.cocons.check.model.list <- function(model.list, data) {
  
  stopifnot("model.list not a list" = is.list(model.list))
  
  if(is.null(model.list$std.dev) || is.null(model.list$scale)){
    stop("scale and std.dev must be specified.")
  }

  if (any(!names(model.list) %in% getOption("cocons.Dictionary"))) {
    stop("aspect names do not match reference ones. Please check getOption(\"cocons.Dictionary\")")
  }
  
  lapply(model.list, FUN = function(x) {
    if (any(!(all.vars(x)) %in% colnames(data))) {
      stop("variable names in model.list do not match data variable names.")
    }
  })
  
  return(0)
}

.cocons.check.info <- function(type, info, model.list, data){
  
  if (is.null(info$smooth.limits) & is.formula(model.list[6]$smooth)) {
    stop("smooth limits not specified.")
  }
  
  if (!is.null(info$lambda.reg)) {
    if(info$lambda.reg < 0){
      stop("lambda.reg must be non-negative.")
    }
  }
  
  if (!is.null(info$lambda.Sigma)) {
    if(info$lambda.Sigma < 0){
      stop("lambda.Sigma must be non-negative.")
    }
  }
  
  if (!is.null(info$lambda.betas)) {
    if(info$lambda.betas < 0){
      stop("lambda.betas must be non-negative.")
    }
  }
    
  if(!is.null(info$smooth.limits)){
    if (info$smooth.limits[1] < 0) {
      stop("lower bound smooth_limit is <= 0 . Should be > 0")
    }
    
    if (info$smooth.limits[1] > info$smooth.limits[2]) {
      stop("lower bound smooth_limit is > upper bound . Should be the opposite.")
    }
  }
  
  if(!is.null(info$smooth.limits)){
    if(is.formula(model.list$smooth) && (info$smooth.limits[1] == info$smooth.limits[2])){
      stop("cannot estimate the smoothness when smooth.limits[1] = smooth.limits[2]")
    }
  }
  
  # sparse
  if (type == "sparse") {
    # taper function
    if (is.null(info$taper) && identical(info$taper, function() {})) {
      stop("taper must be one compact supported function from package spam")
    }
    
    # delta
    if (is.null(info$delta)) {
      stop("taper type requires specifying a delta > 0")
    }
    
    # delta
    if (info$delta < 0) {
      stop("taper argument must be non-negative")
    }
  }
  
  if (type == "dense" && !is.null(info$taper) && !identical(info$taper, function() {})) {
    stop("if type is dense taper should not be specified")
  }
  if (type == "dense" && !is.null(info$delta) && !is.na(info$delta)) {
    stop("if type is dense do not specify delta")
  }
  
  if (!is.null(info$skip.scale)) {
    
    if(!all(info$skip.scale > 0 & info$skip.scale <= dim(data)[2])){
      stop("check skip.scale.")
      
    }
    
  }
  
}

.cocons.check.output <- function(output){
  # output
  if (!(identical(output, list()))) {
    warning("providing an output object. Check optim.coco() to match output list type.")
  }
}

.cocons.check.Hess <- function(Hess) {
  return(0)
}

.cocons.check.ncores <- function(ncores){
  #if(!is.integer(ncores)){stop("ncores needs to be an integer.")}
  if(ncores > parallel::detectCores()){stop("ncores must be less than available cores.")}
  if(ncores < 1){stop("ncores must be >=1")}
}

.cocons.check.boundaries <- function(coco.object, boundaries){
  
  tmp_boundaries <- cocons::getBoundaries(
    x = coco.object, 
    lower.value = -2,
    upper.value = 2)
  
  if(length(tmp_boundaries) != length(boundaries)){
    stop("check boundaries.")
  }
  
  if(any(unlist(lapply(tmp_boundaries,length)) != unlist(lapply(boundaries, length)))){
    stop("check boundaries.")
  }
  
  if(any(unlist(lapply(tmp_boundaries,is.na)))){
    print("NAs in the boundaries not allowed.")
  }

}

.cocons.check.newdataset <- function(newdataset, coco.object){
  
  if(!is.data.frame(newdataset)){
    stop("newdataset must be a data.frame object.")
  }
  
  colnames_nd <- colnames(getDesignMatrix(model.list = coco.object@model.list, data = newdataset)$model.matrix)
  colnames_co <- colnames(getDesignMatrix(model.list = coco.object@model.list, data = coco.object@data)$model.matrix)
  
  if(any(colnames_nd != colnames_co)){
    stop("check colnames of newdataset.")
  }
  
  if(any(is.na(newdataset))){
    stop("NAs in newdataset not allowed.")
  }
  
}

.cocons.check.newlocs <- function(newlocs){
  
  if(!is.matrix(newlocs)){
    stop("newlocs is not a matrix.")
  }
  
  if(dim(newlocs)[2] != 2){
    stop("check dimension of newlocs.")
  }
  
  if(any(is.na(newlocs))){
    stop("newlocs with NAs.")
  }
  
}

.cocons.check.object <- function(object){
  
}

.cocons.check.convergence <- function(output, boundaries){
  
  if(any(boundaries$theta_upper == output$par, na.rm = T) ||
      any(boundaries$theta_lower == output$par, na.rm = T)) {
    warning("at least one of the estimates at the 
              boundaries.")
  }
  
  if(any(output$loginfo[,which(colnames(output$loginfo) == 'fn')] == 1e6)){
    which_ones <- which(output$loginfo[,which(colnames(output$loginfo) == 'fn')] == 1e6)
    warning("ill-posed covariance matrix at evaluation/s ", paste0(which_ones,collapse = ","))
  }
  
  if(output$convergence != 0){
    return(print(output$message))
  }
  
}

.cocons.setDesignMatrixCat <- function(coco.object, designMatrix){
  
  to_not_std <- colnames(coco.object@data)[coco.object@info$skip.scale]
  to_avoid_std <- colnames(designMatrix$model.matrix) %in% to_not_std
  
  tmp_values <- cocons::getScale(designMatrix$model.matrix[,!to_avoid_std])
  
  empty_matrix <- matrix(0, ncol = dim(designMatrix$model.matrix)[2], 
                         nrow = dim(designMatrix$model.matrix)[1])
  
  empty_matrix[, !to_avoid_std] <- tmp_values$std.covs
  empty_matrix[, to_avoid_std] <- designMatrix$model.matrix[, to_avoid_std]
  
  mean_vector_empty <- rep(0,dim(designMatrix$model.matrix)[2])
  mean_sd_empty <- rep(1,dim(designMatrix$model.matrix)[2])
  
  mean_vector_empty[!to_avoid_std] <- tmp_values$mean.vector
  mean_sd_empty[!to_avoid_std] <- tmp_values$sd.vector
  
  tmp_values$mean.vector <- mean_vector_empty
  tmp_values$sd.vector <- mean_sd_empty
  
  return(list('mod_DM' = empty_matrix,
              'tmp_values' = tmp_values))
}

.cocons.getPen <- function(n, lambda, theta_list, smooth.limits){
  
  if(T){
    summ <- lambda[3] * exp(theta_list$scale[1]) * 
      sqrt(((smooth.limits[2]-smooth.limits[1])/ 
              (1 + exp(-theta_list$smooth[1])) + 
              smooth.limits[1])) + sumsmoothlone(x = unlist(theta_list[[1]][-1]), lambda[2])    
  }

  for(ii in 2:6){
    
    lambda_i <- sumsmoothlone(x = unlist(theta_list[[ii]][-1]), lambda[1])
    summ <- summ + lambda_i
    
  }
  
  return(2 * n * summ)
  
}

.cocons.getDelta <- function(n, sigma, topDelta = 1e-9){
  
  return(sigma * topDelta * (1/(1+exp(- (n - 5000)/1000))))
  
}

.cocons.updateNames <- function(namesToUpdate,DesignMatrix){
  
  for(ii in 1:length(namesToUpdate)){
    
    tmp_colnames <- colnames(DesignMatrix)
    
    namesToUpdate[ii] <- paste(sub("[0-9]+$", "", namesToUpdate[ii]), tmp_colnames[as.numeric(sub("^[a-zA-Z\\.]+", "", namesToUpdate[ii]))] )
    
  }
  
  return(namesToUpdate)
  
}


.cocons.update.coco.first.step <- function(coco.object, output, boundaries){
  
  # Update model structure
  
  par_pos <- getDesignMatrix(coco.object@model.list,data = coco.object@data)$par.pos
  
  number_parameters_models <- lapply(par_pos,sum)
  
  new_formulas_list <- coco.object@model.list
  
  coco.object@output <- output
  
  parss <- getEstims(coco.object)
  
  # What if one estimate is exactly 0 ?
  # Check models formulas keeping environment
  
  for(ii in 1:7){
    
    if(!is.formula(coco.object@model.list[[ii]])){
      next
    }
    
    to_zero <- which(abs(parss[[ii]][par_pos[[ii]]]) <= coco.object@info$sparse.point)
    
    if(length(to_zero) == 0){
      next
    }
    
    to_zero_covs <- sum(abs(parss[[ii]][par_pos[[ii]]]) <= coco.object@info$sparse.point)
    
    if( to_zero_covs == (number_parameters_models[[ii]] - 1) || (to_zero_covs == (number_parameters_models[[ii]]))){
      new_formulas_list[[ii]] <- stats::as.formula("~1", env = globalenv())
      next
    }
    
    # Update formula
    
    # If the only parameter to drop is the intercept
    if(length(to_zero) == 1 && to_zero == 1){
      next
    }
    
    new_terms <- stats::drop.terms(stats::terms(coco.object@model.list[[ii]]), dropx = (to_zero - 1), keep.response = TRUE)
    
    if(attr(new_terms,"intercept") & (length(attr(new_terms,"term.labels")) > 0)){
      new_terms <- c("1",attr(new_terms, "term.labels"))
    } else{
      new_terms <- attr(new_terms, "term.labels")
    }
    
    new_formula <- stats::reformulate(new_terms, response = NULL, env = globalenv())
    new_formulas_list[[ii]] <- new_formula
  
  }
  
  index_where <- 0
  
  vector_logical <- logical(length = length(output$par))
  
  for(ii in 1:7){
    
    if(!is.formula(coco.object@model.list[[ii]])){
      next
    }
    
    to_zero <- which(abs(parss[[ii]][par_pos[[ii]]]) <= coco.object@info$sparse.point)
    
    if(any(to_zero == 1) && (length(to_zero) > 1)){
      to_zero <- to_zero[-1]
    }
    
    if(length(to_zero) == 1 && to_zero == 1){
      index_where <- index_where + number_parameters_models[[ii]]
      next
    }
    
    index_to_place <- (index_where + 1 ) : (index_where + number_parameters_models[[ii]])
    
    vector_logical[index_to_place][to_zero] <- TRUE
    
    index_where <- index_where + number_parameters_models[[ii]]
    
  }
  
  boundaries_new <- lapply(boundaries, function(x) x[!vector_logical])
  
  coco.object@info$boundaries <- boundaries_new
  
  coco.object@model.list <- new_formulas_list
  
  return(coco.object)
}
