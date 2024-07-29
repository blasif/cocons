#' coco Class
#' @description Creates a coco S4 object.
#' @details 
#' This S4 class is the centerpiece of the cocons package. Two types of coco objects can be created. One is a "dense" type, meaning that the associated nonstationary covariance functions are dense, or 
#' a "sparse" object, which, in combination with the Tapering approach, induce zeroes in the covariance matrix to make it sparse and to unlock a set of efficient algorithms to speed up estimation and prediction routines.
#' 
#' Another important component of the coco S4 class is the model.list specification, involving different formulas provided as a list, where each of them specifies a source of nonstationarity, based on covariates. It involves "trend" for the spatial trend,
#' the "std. dev" for the marginal standard deviation, "scale", "aniso" and "tilt", each of them shaping specific aspects of the local spatial geometrically anisotropy structure,
#' "smooth" handling local smoothness, and "nugget" handling the local nugget effect.
#' 
#' Lastly, arguments for the info list argument involve: 
#' - 'lambda': a positive scalar specifying the regularization parameter.
#' - 'smooth_limits': specifying the allowed range of variation for the spatially varying smoothness.
#' - 'taper': specifying the desired taper function from the spam package (only for 'sparse' coco objects).
#' - 'delta': specifying the taper range/scale (only for 'sparse' coco objects).
#' - 'cat.vars': index of those variables in the data object which are categorical or should not be scaled during the optimization.
#' 
#' @usage coco(type, data, locs, z, model.list, info, output = list())
#' @param type One of two available types "dense" or "sparse". See description.
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
      info$smooth_limits <- c(model.list$smooth[1],
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
