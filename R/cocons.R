#' coco Class
#' @description Creates a coco S4 object.
#' @details 
#' This S4 class is the centerpiece of the \pkg{cocons} package. Two types of coco objects can be created. One is a \code{'dense'} type, meaning that the associated nonstationary covariance functions are dense, or 
#' a \code{'sparse'} object, which, in combination with the Tapering approach, induce zeroes in the covariance matrix to make it sparse and to unlock a set of efficient algorithms to speed up estimation and prediction routines.
#' 
#' An important component of the coco S4 class is the \code{model.list} specification, involving different formulas provided as a list, where each of them specifies a covariate-based parametric model for a specific source of nonstationarity. It involves \code{"trend"} for the spatial trend,
#' the \code{"std. dev"} for the marginal standard deviation, \code{"scale"}, \code{"aniso"} and \code{"tilt"}, each of them shaping specific aspects of the local spatial geometrically anisotropy structure,
#' \code{"smooth"} handling local smoothness, and \code{"nugget"} handling the local nugget effect.
#' 
#' Lastly, arguments for the info list argument involve: 
#' - \code{'lambda'}: a positive scalar specifying the regularization parameter.
#' - \code{'smooth.limits'}: specifying the allowed range of variation for the spatially varying smoothness.
#' - \code{'taper'}: specifying the desired taper function from the spam package (only for 'sparse' coco objects).
#' - \code{'delta'}: specifying the taper range/scale (only for 'sparse' coco objects).
#' - \code{'cat.vars'}: index of those variables in the data object which are categorical or should not be scaled during the optimization.
#' 
#' @usage coco(type, data, locs, z, model.list, info, output = list())
#' @param type One of two available types \code{'dense'} or \code{'sparse'}. See description.
#' @param data A data.frame with covariates information, where colnames(data) matches model.list specification.
#' @param locs A matrix with locs matching data.
#' @param z A matrix of n x r response realizations, one realization per column. When considering only one realization, a vector can also be provided.
#' @param model.list A list specifying a model for each aspect of the spatial structure.
#' @param info A list specifying characteristics of the coco object.
#' @param output Empty or an output from optimparallel output, including as well boundaries.
#' 
#' @returns Creates a coco object.
#' 
#' @author Federico Blasi
#' @seealso [spam::cov.wend1()]
coco <- function(type,
                  data,
                  locs,
                  z,
                  model.list,
                  info = list(),
                  output = list()) {
  
  .cocons.check.type(type)
  .cocons.check.data(data)
  .cocons.check.locs(locs)
  
  if(is.vector(z)){
    z <- matrix(z, ncol = 1)
  }
  
  .cocons.check.z(z)
  
  .cocons.check.model.list(model.list,
                          data)
  .cocons.check.info(type = type, 
                    info = info,
                    model.list = model.list)
  
  .cocons.check.output(output)
  
  if(is.null(info$lambda)){
    info$lambda <- 0
  }
  
  if(T){
    if(!is.formula(model.list$smooth)){
      info$smooth.limits <- c(model.list$smooth[1],
                              model.list$smooth[1])
    }
  }
  
  methods::new(
    Class = "coco",
    type = type,
    data = data,
    locs = locs,
    z = z,
    model.list = model.list,
    info = info,
    output = output
  )
  
}
