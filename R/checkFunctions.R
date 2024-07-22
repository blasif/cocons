
# Internal functions

# ADD checks for NA's!!!!
# ADD checks for different sizes of locs, data, and z !!!!!!

.coco.update.optim.control <- function(optim.control){
  
  to_update <- getOption("coco.Optim.Control")
  
  if(is.null(optim.control$parallel$forward)){
    optim.control$parallel$forward <- to_update$parallel$forward
  }
  
  if(is.null(optim.control$parallel$loginfo)){
    optim.control$parallel$loginfo <- to_update$parallel$loginfo
  }
  
  if(is.null(optim.control$control$factr)){
    optim.control$control$factr <- to_update$control$factr
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
  
  optim.control$control$hessian <- FALSE

  return(optim.control)
  
}

.coco.DrawEllipsoid <- function(alpha_i, r, rho, loc, factr) {
  
  createEllipse <- function(center, a, b, angle = 0, steps = 100, factr = 0.1) {
    theta <- seq(0, 2 * pi, length.out = steps)
    ellipse_points <- cbind(a * cos(theta), b * sin(theta))
    ellipse_points <- ellipse_points %*% matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
    ellipse_points[, 1] <- factr * ellipse_points[, 1] + center[1]
    ellipse_points[, 2] <- factr * ellipse_points[, 2] + center[2]
    
    return(ellipse_points)
  }
  
  if (alpha_i == pi / 2) {
    print("linear independent eigenvectors")
    A_stuff <- sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2)
    
    e_1 <- ((r + 1) + sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
    e_2 <- ((r + 1) - sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
    
    magnitude_vector_one <- sqrt((0)^2 + 1^2)
    y_length <- sqrt(max(e_2, e_1)) / magnitude_vector_one
    
    graphics::arrows(
      x0 = 0, x1 = ifelse(r > 1, 0, y_length),
      y0 = 0, y1 = ifelse(r > 1, y_length, 0), lwd = 2, lty = 1,
      cex = 0.5, angle = 5, length = 0.1
    )
    
    magnitude_vector_two <- sqrt((1)^2 + (0)^2)
    
    y_length <- sqrt(min(e_1, e_2)) / magnitude_vector_two
    
    graphics::arrows(
      x0 = 0, x1 = ifelse(r > 1, y_length, 0),
      y0 = 0, y1 = ifelse(r > 1, 0, y_length), lwd = 2, lty = 1,
      cex = 0.5, angle = 5, length = 0.1
    )
    
    ellipse_points <- createEllipse(loc,
                                    a = sqrt(ifelse(r > 1, (min(e_1, e_2)), max(e_1, e_2))),
                                    b = sqrt(ifelse(r > 1, (max(e_1, e_2)), min(e_1, e_2))), 0, steps = 2000
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
  
  graphics::arrows(
    x0 = loc[1], x1 = loc[1] + y_length * a_a_1,
    y0 = loc[2], y1 = loc[2] + y_length * a_b_1, lwd = 2, lty = 1,
    cex = 0.5, angle = 15, length = 0.1, col = "black"
  )
  
  magnitude_vector_two <- sqrt((e_b_1)^2 + (e_b_2)^2)
  
  y_length <- factr * sqrt(e_1) / magnitude_vector_two
  
  graphics::arrows(
    x0 = loc[1], x1 = loc[1] + y_length * e_b_1,
    y0 = loc[2], y1 = loc[2] + y_length * e_b_2, lwd = 2, lty = 1,
    cex = 0.5, angle = 15, length = 0.1, col = "black"
  )
  
  ellipse_points <- createEllipse(
    center = loc, a = sqrt((max(e_1, e_2))),
    b = sqrt((min(e_1, e_2))), atan2(x = e_b_1, y = e_b_2),
    steps = 2000, factr = factr
  )
  
  graphics::lines(ellipse_points, lty = 3, col = "black")
}

.coco.check.type_pred <- function(type){}

.coco.check.type <- function(type) {
  if (!(type %in% c("sparse", "dense"))) {
    stop("check type")
  }
  return(0)
}

.coco.check.data <- function(data) {
  if (!is.data.frame(data)) {
    stop("data must be provided as data.frame")
  }
  
  if (is.null(colnames(data))) {
    stop("data does not has colnames, which are needed for building `par.pos` ")
  }
  return(0)
}

.coco.check.locs <- function(locs) {
  if (!is.matrix(locs)) {
    stop("locs should be provided as matrix")
  }
  return(0)
}

.coco.check.z <- function(z) {
  if (is.null(z)) {
    warning("z not provided. Expecting to simulate with this coco object.")
  } else {
    if ( !is.matrix(z)) {
      stop("z is not a matrix.")
    }
  }
  return(0)
}

.coco.check.model.list <- function(model.list, data) {
  if (!is.list(model.list)) {
    stop("model.list not a list")
  }
  
  if (any(names(model.list) != getOption("coco.Dictionary"))) {
    stop("aspect names do not match reference ones. Please check getOption(\"coco.Dictionary\")")
  }
  
  lapply(model.list, FUN = function(x) {
    if (any(!(all.vars(x)) %in% colnames(data))) {
      stop("variable names in list_formula do not match data variable names.")
    }
  })
  
  return(0)
}

# ADDED model.list here to check
.coco.check.info <- function(type, info, model.list){
  
  if (is.null(info$smooth_limits)) {
    warning("smooth limits not specified. Maybe considering fixed smoothness?")
  }
  
  if (!is.null(info$lambda)) {
    if(info$lambda < 0){
      stop("lambda must be non-negative.")
    }
  }
    
  if(!is.null(info$smooth_limits)){
    if (info$smooth_limits[1] < 0) {
      stop("lower bound smooth_limit is < 0 . Should be > 0")
    }
    
    if (info$smooth_limits[1] > info$smooth_limits[2]) {
      stop("lower bound smooth_limit is > upper bound . Should be the opposite.")
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
  
}

.coco.check.output <- function(output){
  # output
  if (!(identical(output, list()))) {
    warning("providing an output object. Check optim.coco() to match output list type.")
  }
}

.coco.check.Hess <- function(Hess) {
  return(0)
}

.coco.check.ncores <- function(ncores){
  #if(!is.integer(ncores)){stop('ncores needs to be an integer.')}
  if(ncores > parallel::detectCores()){stop('ncores must be less than available cores.')}
  if(ncores < 1){stop('ncores must be >=1')}
}

.coco.check.boundaries <- function(boundaries){
  
}

.coco.check.newdataset <- function(newdataset){
  
}

.coco.check.newlocs <- function(newlocs){
  
}

.coco.check.object <- function(object){
  
}
