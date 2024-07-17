#' svcov Class
#' @description Creates a svcov object
#' @usage svcov(type, model.list, data, locs, z, 
#' taper = function(){}, delta = NA_real_, output = list())
#' @slot type One of two available types "dense" or "sparse". See description.
#' @slot model.list A list specyfing a model for each aspect of the spatial structure.
#' @slot data A data.frame with covariates information, where colnames(data) matches model.list specification
#' @slot locs a matrix with locs matching data
#' @slot z A vector with response values
#' @slot taper If "sparse", a compact supported taper function from spam package
#' @slot delta a taper range
#' @slot output empty or an output from optimparallel output, including as well boundaries 
#' @returns Creates an svcov object
#' @author Federico Blasi
#' @seealso [spam::cov.wend1()]
#' 
svcov <- function(type,
                  data,
                  locs,
                  z,
                  model.list,
                  info = list(),
                  output = list()) {
  
  .svcov.check.type(type)
  .svcov.check.data(data)
  .svcov.check.locs(locs)
  
  if(is.vector(z)){
    z <- matrix(z, ncol = 1)
  }
  
  .svcov.check.z(z)
  
  .svcov.check.model.list(model.list,
                          data)
  .svcov.check.info(type = type, 
                    info = info,
                    model.list = model.list)
  
  .svcov.check.output(output)
  
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
    Class = "svcov",
    type = type,
    data = data,
    locs = locs,
    z = z,
    model.list = model.list,
    info = info,
    output = output
  )
  
}
